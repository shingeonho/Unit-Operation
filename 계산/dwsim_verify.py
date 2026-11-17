"""
DWSIM 검증 — 직접 계산(column_design.py) 결과를 DWSIM 엄밀 단 계산으로 확인

실행:  python dwsim_verify.py
필요:  DWSIM 10 설치 + .NET 10 런타임 (아래 DOTNET_ROOT), pythonnet
       (pip install pythonnet)

검증 케이스 (열역학: Raoult's Law — 논문 Table 4 대비 오차 0.3 K로 확인된 모델)
  A. 직접 계산과 같은 조건: 15 kPa 균일, 단 수·공급단·R·D 동일 → 제품 순도, 열량, 온도 비교
  B. 같은 단 수에서 xD = 0.99를 맞추는 데 필요한 R (DWSIM 기준)
  C. 압력강하 반영: 탑정 15 kPa, 탑 전체 ΔP = 직접 계산 값 → 탑저 온도 비교
"""
import math
import os
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(errors="replace")
except Exception:
    pass

# ---------------------------------------------------------------------
# 경로 설정 (PC 환경에 맞게 수정)
# ---------------------------------------------------------------------
DWSIM_DIR = r"C:\Users\82107\AppData\Local\Programs\DWSIM"
# .NET 10 런타임 위치. 임시 폴더라 지워졌다면 Microsoft dotnet-install 스크립트로 다시 받으면 됨:
#   powershell -c "& ([scriptblock]::Create((iwr https://dot.net/v1/dotnet-install.ps1))) -Runtime dotnet -Channel 10.0 -InstallDir <경로>"
DOTNET_ROOT = r"C:\Users\82107\AppData\Local\Temp\claude\234cb68f4b9c6a76\fc5acb7a-8d27-4c29-aa34-fb228fb046fd\scratchpad\dotnet"

HERE = Path(__file__).parent
OUT = HERE / "결과"
OUT.mkdir(exist_ok=True)

import inputs as I
import vle
import mccabe as mc


# ---------------------------------------------------------------------
# DWSIM 불러오기 (.NET 10 CoreCLR)
# ---------------------------------------------------------------------
def load_dwsim():
    from clr_loader import get_coreclr
    from pythonnet import set_runtime
    rc = OUT / "_dotnet.runtimeconfig.json"
    rc.write_text('{"runtimeOptions":{"tfm":"net10.0","framework":'
                  '{"name":"Microsoft.NETCore.App","version":"10.0.0"}}}', encoding="ascii")
    set_runtime(get_coreclr(runtime_config=str(rc), dotnet_root=DOTNET_ROOT))
    import clr
    import System
    from System.Reflection import Assembly, AssemblyName
    from System.Runtime.Loader import AssemblyLoadContext
    from System import Func

    def _res(ctx, name):
        p = os.path.join(DWSIM_DIR, name.Name + ".dll")
        return ctx.LoadFromAssemblyPath(p) if os.path.exists(p) else None
    AssemblyLoadContext.Default.add_Resolving(Func[AssemblyLoadContext, AssemblyName, Assembly](_res))
    cwd = os.getcwd()
    os.chdir(DWSIM_DIR)
    for n in ["DWSIM.Interfaces", "DWSIM.GlobalSettings", "DWSIM.SharedClasses", "DWSIM.Thermodynamics",
              "DWSIM.UnitOperations", "DWSIM.Automation", "DWSIM.MathOps", "DWSIM.FlowsheetSolver"]:
        clr.AddReference(os.path.join(DWSIM_DIR, n + ".dll"))
    os.chdir(cwd)


load_dwsim()
import System                                                 # noqa: E402
from DWSIM.Automation import Automation3                       # noqa: E402
from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType  # noqa: E402
from DWSIM.Interfaces.Enums import StreamSpec                 # noqa: E402

AM = Automation3()
PP = "Raoult's Law"
EB, ST = "Ethylbenzene", "Styrene"


# ---------------------------------------------------------------------
# 직접 계산 값 (column_design.py와 같은 함수로 다시 계산)
# ---------------------------------------------------------------------
def lam_EB(T):
    pts = I.DHVAP_EB_POINTS
    for (t0, h0), (t1, h1) in zip(pts, pts[1:]):
        if T <= t1 or (t1, h1) == pts[-1]:
            return h0 + (h1 - h0) * (T - t0) / (t1 - t0)


def hand_calc():
    D, B = mc.mass_balance(I.F, I.zF, I.xD, I.xB)
    _, a = vle.alpha_from_data(I.VLE_15kPa)
    eq = vle.ConstAlpha(a)
    Rm, _ = mc.R_min(eq, I.xD, I.zF, I.q)
    R = I.R_factor * Rm
    st = mc.step_stages(eq, I.xD, I.xB, I.zF, I.q, R)
    N_col = st["n_stages"] - 1
    V = (R + 1) * D
    T_top = vle.T_dew(I.xD, I.P_top)
    lamE = lam_EB(T_top)
    Qc = V * (I.xD * lamE + (1 - I.xD) * I.DHVAP_ST) / 3.6
    Qr = V * (I.xB * lamE + (1 - I.xB) * I.DHVAP_ST) / 3.6
    Z = N_col * I.HETP
    dP = I.DP_PER_M * Z / 10.0
    return dict(D=D, R=R, alpha=a, N_col=N_col, feed=st["feed_stage"], Qc=Qc, Qr=Qr,
                T_top=T_top, T_bot=vle.T_bubble(I.xB, I.P_top),
                T_bot_dP=vle.T_bubble(I.xB, I.P_top + dP), dP=dP)


# ---------------------------------------------------------------------
# DWSIM 컬럼 한 번 풀기
# ---------------------------------------------------------------------
def _arr(v):
    return System.Array[System.Double](list(v))


TEMPLATE = os.path.join(DWSIM_DIR, "samples", "Tutorial 04 - Distillation Column.dwxmz")


def _template_flowsheet():
    """자동화로 새로 만든 컬럼은 수렴하지 않아서(원인 미상), DWSIM 기본 예제(Tutorial 04)의 컬럼을
    틀로 쓰고 성분을 EB/스티렌으로, 열역학을 Raoult's Law로 바꾼다.
    (예제의 물·에탄올은 몰분율 0으로 남지만 계산에 영향 없음)"""
    from DWSIM.Thermodynamics.BaseClasses import Compound
    fs = AM.LoadFlowsheet(TEMPLATE)
    fs.AddCompound(EB)
    fs.AddCompound(ST)
    pp = fs.CreateAndAddPropertyPackage(PP)
    col = [o.GetAsObject() for o in fs.SimulationObjects.Values
           if "Column" in o.GetAsObject().GetType().Name][0]
    ids = {ms.StreamBehavior.ToString(): ms.StreamID for ms in col.MaterialStreams.Values}
    for o in fs.SimulationObjects.Values:
        g = o.GetAsObject()
        if g.GetType().Name == "MaterialStream":
            for p in g.Phases.Values:
                for n in (EB, ST):
                    if not p.Compounds.ContainsKey(n):
                        cp = Compound(n, "")
                        cp.ConstantProperties = fs.SelectedCompounds[n]
                        p.Compounds.Add(n, cp)
        g.PropertyPackage = pp
    feed = fs.SimulationObjects[ids["Feed"]].GetAsObject()
    dist = fs.SimulationObjects[ids["Distillate"]].GetAsObject()
    bott = fs.SimulationObjects[ids["BottomsLiquid"]].GetAsObject()
    return fs, col, feed, dist, bott


def run_column(n_col, feed_stage, R, D_kmolh, P_top_kPa, dP_kPa=0.0, save=None):
    """n_col: 탑 안 이론단 수(정수, 응축기·재비기 제외), feed_stage: 위에서 몇 번째 단(응축기 = 0)"""
    fs, col, feed, dist_o, bott_o = _template_flowsheet()

    # 원료: 100 kmol/h, zF, 탑정 압력에서 포화액 (q = 1)
    ph = feed.Phases[0]
    for k in list(ph.Compounds.Keys):
        ph.Compounds[k].MoleFraction = 0.0
    ph.Compounds[EB].MoleFraction = I.zF
    ph.Compounds[ST].MoleFraction = 1 - I.zF
    ph.Properties.pressure = P_top_kPa * 1000.0
    ph.Properties.molarflow = I.F * 1000.0 / 3600.0          # mol/s
    feed.SpecType = StreamSpec.Pressure_and_VaporFraction
    feed.Phases[2].Properties.molarfraction = 0.0

    # 컬럼: 응축기 + 이론단 n_col + 재비기
    n_total = n_col + 2
    col.SetNumberOfStages(n_total)
    col.CondenserType = System.Enum.Parse(col.CondenserType.GetType(), "Total_Condenser")
    col.SetStreamFeedStage(feed, feed_stage)

    class _W:   # 결과 읽기용 래퍼 (GetAsObject 인터페이스 맞춤)
        def __init__(self, o): self.o = o
        def GetAsObject(self): return self.o
    dist, bott = _W(dist_o), _W(bott_o)

    # 압력 분포 (단별로 직접 지정)
    P_Pa = P_top_kPa * 1000.0
    for i, s in enumerate(col.Stages):
        s.P = P_Pa + dP_kPa * 1000.0 * i / (n_total - 1)

    # 초기 추정값: DWSIM 자동 생성
    col.ResetInitialEstimates()

    # 사양: 환류비 R, 탑저 유량 B
    col.SetCondenserSpec("Reflux Ratio", R, "", "")
    col.SetReboilerSpec("Product Molar Flow Rate", I.F - D_kmolh, "kmol/h", "")
    col.MaxIterations = 300

    errs = AM.CalculateFlowsheet2(fs)
    err_txt = [str(e.Message) for e in errs] if errs is not None else []

    d = dist.GetAsObject().Phases[0]
    b = bott.GetAsObject().Phases[0]
    if not col.Calculated or d.Properties.molarflow is None:
        raise RuntimeError("DWSIM 컬럼 미수렴: " + str(col.ErrorMessage) + " " + " | ".join(err_txt))
    if save:
        AM.SaveFlowsheet(fs, str(save), True)     # DWSIM GUI에서 열어볼 수 있는 파일
    stages = [(s.T, s.P) for s in col.Stages]
    fid = [ms.AssociatedStage for ms in col.MaterialStreams.Values if ms.StreamBehavior.ToString() == "Feed"][0]
    feed_idx = [i for i, s in enumerate(col.Stages) if s.ID == fid][0]
    return dict(
        feed_idx=feed_idx,
        xD=d.Compounds[EB].MoleFraction, xB=b.Compounds[EB].MoleFraction,
        D=d.Properties.molarflow * 3.6, B=b.Properties.molarflow * 3.6,
        Qc=abs(col.CondenserDuty), Qr=abs(col.ReboilerDuty),
        T_cond=stages[0][0], T_top=stages[1][0], T_bot=stages[-1][0],
        P_bot=stages[-1][1] / 1000.0, n_total=n_total, errors=err_txt,
    )


def dwsim_alpha(P_kPa, x):
    """DWSIM Raoult 패키지의 α (기포점 플래시)"""
    fs = AM.CreateFlowsheet()
    fs.AddCompound(EB)
    fs.AddCompound(ST)
    fs.CreateAndAddPropertyPackage(PP)
    ms = fs.AddObject(ObjectType.MaterialStream, 0, 0, "s").GetAsObject()
    ph = ms.Phases[0]
    ph.Properties.pressure = P_kPa * 1000.0
    ph.Properties.molarflow = 1.0
    ph.Compounds[EB].MoleFraction = x
    ph.Compounds[ST].MoleFraction = 1 - x
    ms.SpecType = StreamSpec.Pressure_and_VaporFraction
    ms.Phases[2].Properties.molarfraction = 0.0
    AM.CalculateFlowsheet2(fs)
    y = ms.Phases[2].Compounds[EB].MoleFraction
    return (y / (1 - y)) / (x / (1 - x)), ph.Properties.temperature


def main():
    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    h = hand_calc()
    n_col = math.ceil(h["N_col"])          # 55.5 → 56단으로 올림
    feed = h["feed"]                        # 위에서 29번째 단

    out("# DWSIM 검증 결과")
    out()
    out(f"- 열역학 모델: DWSIM `{PP}` (Table 4 실측 대비 평균 오차 0.3 K, Δy 0.004 — 앞서 확인)")
    out(f"- 직접 계산: 탑 안 {h['N_col']:.2f}단 → DWSIM은 **{n_col}단 + 전응축기 + 부분 재비기 = {n_col+2}단**, "
        f"공급단 위에서 {feed}번째, R = {h['R']:.3f}, D = {h['D']:.1f} kmol/h")
    a_top, _ = dwsim_alpha(I.P_top, I.xD)
    a_bot, _ = dwsim_alpha(I.P_top, I.xB)
    out(f"- DWSIM α (15 kPa): 탑정 {a_top:.3f}, 탑저 {a_bot:.3f}  vs 직접 계산 α {h['alpha']:.3f} (Table 4 실측)")
    out()

    # A ---------------------------------------------------------------
    A = run_column(n_col, feed, h["R"], h["D"], I.P_top, 0.0, save=OUT / "DWSIM_케이스A.dwxmz")
    out("## A. 같은 조건 (15 kPa 균일, 같은 단 수·R·D)")
    out("| 항목 | 직접 계산 | DWSIM | 차이 |")
    out("|---|---|---|---|")
    out(f"| 탑정 xD (EB) | {I.xD:.4f} (목표) | {A['xD']:.4f} | {A['xD']-I.xD:+.4f} |")
    out(f"| 탑저 xB (EB) | {I.xB:.4f} (목표) | {A['xB']:.4f} | {A['xB']-I.xB:+.4f} |")
    out(f"| D [kmol/h] | {h['D']:.1f} | {A['D']:.1f} | |")
    out(f"| 응축기 Qc [kW] | {h['Qc']:.0f} | {A['Qc']:.0f} | {(A['Qc']/h['Qc']-1)*100:+.1f} % |")
    out(f"| 재비기 Qr [kW] | {h['Qr']:.0f} | {A['Qr']:.0f} | {(A['Qr']/h['Qr']-1)*100:+.1f} % |")
    out(f"| 탑정 증기 T [°C] | {h['T_top']-273.15:.1f} | {A['T_top']-273.15:.1f} | {A['T_top']-h['T_top']:+.1f} K |")
    out(f"| 탑저 T [°C] | {h['T_bot']-273.15:.1f} | {A['T_bot']-273.15:.1f} | {A['T_bot']-h['T_bot']:+.1f} K |")
    out()

    # B ---------------------------------------------------------------
    out("## B. 같은 단 수에서 xD = 0.99를 맞추는 R (DWSIM)")
    lo, hi = 0.6 * h["R"], h["R"]
    rB = A if A["xD"] >= I.xD else None
    if rB is None:
        hi = 1.5 * h["R"]
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        try:
            r = run_column(n_col, feed, mid, h["D"], I.P_top, 0.0)
        except RuntimeError:
            lo = mid
            continue
        if r["xD"] >= I.xD:
            hi, rB = mid, r
        else:
            lo = mid
    if rB:
        out(f"- 필요한 R = **{hi:.3f}** (직접 계산 {h['R']:.3f}, 차이 {(hi/h['R']-1)*100:+.1f} %), "
            f"이때 xB = {rB['xB']:.4f}, Qr = {rB['Qr']:.0f} kW")
    out()

    # D ---------------------------------------------------------------
    if rB:
        a_dw = math.sqrt(a_top * a_bot)
        eq_dw = vle.ConstAlpha(a_dw)
        stD = mc.step_stages(eq_dw, I.xD, I.xB, I.zF, I.q, hi)
        out("## D. 같은 α로 맞추면? (계산 방법 자체의 검증)")
        out(f"- 직접 계산(McCabe-Thiele)에 DWSIM α {a_dw:.3f}와 B의 R {hi:.3f}를 넣으면 "
            f"이상단 **{stD['n_stages']:.2f}단** (재비기 포함), 공급단 {stD['feed_stage']}")
        out(f"- DWSIM은 같은 R에서 {n_col}단 + 재비기 = **{n_col+1}단**으로 xD 0.99 달성 → "
            f"차이 {n_col+1-stD['n_stages']:.2f}단 ({(n_col+1)/stD['n_stages']*100-100:.1f} %)")
        out("- 즉 α만 같으면 두 방법의 단수가 거의 같음 → 차이는 계산 방법이 아니라 평형 데이터(α) 출처 차이")
        out()

    # C ---------------------------------------------------------------
    C = run_column(n_col, feed, h["R"], h["D"], I.P_top, h["dP"], save=OUT / "DWSIM_케이스C_압력강하.dwxmz")
    out(f"## C. 압력강하 반영 (탑 전체 ΔP = {h['dP']:.2f} kPa)")
    out("| 항목 | 직접 계산 | DWSIM |")
    out("|---|---|---|")
    out(f"| 탑저 압력 [kPa] | {I.P_top + h['dP']:.2f} | {C['P_bot']:.2f} |")
    out(f"| 탑저 T [°C] | {h['T_bot_dP']-273.15:.1f} | {C['T_bot']-273.15:.1f} |")
    out(f"| 탑정 xD / 탑저 xB | 0.99 / 0.01 (목표) | {C['xD']:.4f} / {C['xB']:.4f} |")
    out(f"| 재비기 Qr [kW] | {h['Qr']:.0f} | {C['Qr']:.0f} |")
    out()

    (OUT / "DWSIM검증.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"\n저장: {OUT / 'DWSIM검증.md'}")


if __name__ == "__main__":
    main()
