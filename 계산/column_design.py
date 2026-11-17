"""
EB–스티렌 감압 충전 증류탑 설계 (직접 계산, 강의 이론 기반)

실행:  python column_design.py      (VS Code에서 이 파일을 열고 ▶ 실행)
입력:  inputs.py      결과:  결과/ 폴더 (요약 md, 그림 png)

계산 순서 (강의 근거)
  0. 방법 검증: Example 21.2 재현                         (SP6 p30~36)
  1. 물질수지                                              (SP6 p9)
  2. 평형: α = Table 4 실측(15 kPa) 기하평균 [주 계산],
           Raoult α = P_EB'/P_ST' (Table 2) [비교]           (SP6 p4, p25)
  3. Rmin, Nmin(Fenske + 전환류 작도), R = 1.2 Rmin, 이상단 수 (SP6 p21~29)
  4. 유량, 응축기·재비기 열량, 정몰유량 가정 점검           (SP6 p24, SP7 p1~4)
  5. 범람 속도(Fig 18.8), 직경 (탑정·탑저 중 큰 값)          (SP3 p11~14)
  6. 충전 높이 Z = N × HETP (HETP ≈ H_Oy)                  (McCabe 21.6)
  7. 압력강하로 탑저 압력·온도 점검 (Salem 실측 mbar/m)       (설계 후 점검)
"""
import math
import sys
from pathlib import Path
import matplotlib

try:   # 한글 Windows 콘솔(cp949)에서 특수기호 출력 오류 방지
    sys.stdout.reconfigure(errors="replace")
except Exception:
    pass
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import inputs as I
import vle
import mccabe as mc
import example21_2_test

OUT = Path(__file__).parent / "결과"
OUT.mkdir(exist_ok=True)
plt.rcParams["font.family"] = "Malgun Gothic"
plt.rcParams["axes.unicode_minus"] = False

lines = []          # 요약 md에 쓸 줄


def out(s=""):
    print(s)
    lines.append(s)


def lam_EB(T):
    """NIST 열량계 값(Svoboda 1982) 선형 보간 [kJ/mol]"""
    pts = I.DHVAP_EB_POINTS
    for (t0, h0), (t1, h1) in zip(pts, pts[1:]):
        if T <= t1 or (t1, h1) == pts[-1]:
            return h0 + (h1 - h0) * (T - t0) / (t1 - t0)


def interp_fig18_8(X):
    """Fig 18.8 Mellapak 250Y: Y를 log X에 대해 선형 보간 (범위 밖은 끝값 유지)"""
    pts = I.FIG18_8_MELLAPAK250Y
    if X <= pts[0][0]:
        return pts[0][1], "범위 아래(끝값 사용)"
    if X >= pts[-1][0]:
        return pts[-1][1], "범위 위(끝값 사용)"
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= X <= x1:
            t = (math.log10(X) - math.log10(x0)) / (math.log10(x1) - math.log10(x0))
            return y0 + t * (y1 - y0), "판독 범위 내"


def main():
    # ------------------------------------------------------------------
    out("# EB–스티렌 감압 충전탑 직접 계산 결과")
    out()
    t = example21_2_test.run(verbose=True)
    out("## 0. 방법 검증 (Example 21.2)")
    out(f"- D = {t['D']:.1f} kmol/h (교재 153.4), 이상단 {t['n']:.2f} (교재 11, 재비기 포함), "
        f"공급단 {t['feed']} (교재 7) → **{'±1단 이내 일치' if t['ok'] else '확인 필요'}**")
    out("- 차이 원인: 교재는 실측 평형곡선(Fig 21.14), 이 코드는 α = 2.5 일정 근사")
    out()

    # 1. 물질수지 --------------------------------------------------------
    D, B = mc.mass_balance(I.F, I.zF, I.xD, I.xB)
    out("## 1. 설계 조건과 물질수지")
    out(f"- F = {I.F:.1f} kmol/h, zF = {I.zF}, xD = {I.xD}, xB = {I.xB} (EB 몰분율), q = {I.q}")
    out(f"- 탑정 압력 {I.P_top} kPa, R = {I.R_factor} × Rmin, 충전물 {I.packing}, HETP {I.HETP} m, 범람 {I.flood_frac*100:.0f} %")
    out(f"- **D = {D:.2f} kmol/h, B = {B:.2f} kmol/h**")
    out()

    # 2. 평형 ------------------------------------------------------------
    P = I.P_top
    al_pts, a_exp = vle.alpha_from_data(I.VLE_15kPa)
    T_top = vle.T_dew(I.xD, P)          # 탑정 증기(= 전응축기 입구) 이슬점
    T_cond = vle.T_bubble(I.xD, P)      # 환류액(포화액) 온도
    T_bot = vle.T_bubble(I.xB, P)       # 재비기 액체 기포점
    a_top, a_bot = vle.alpha_raoult(T_top), vle.alpha_raoult(T_bot)
    a_raoult_gm = math.sqrt(a_top * a_bot)
    out("## 2. 기액평형 (15 kPa)")
    out(f"- Table 4 실측 α: 범위 {min(al_pts):.3f} ~ {max(al_pts):.3f}, **기하평균 {a_exp:.3f}** ← 주 계산에 사용")
    out(f"- Raoult(Table 2 증기압) α: 탑정 {a_top:.3f} @ {T_top:.2f} K, 탑저 {a_bot:.3f} @ {T_bot:.2f} K, 기하평균 {a_raoult_gm:.3f}")
    # Raoult 모델이 Table 4 실측을 얼마나 재현하는지
    dT = [abs(vle.T_bubble(x, P) - T) for T, x, y in I.VLE_15kPa]
    dy = [abs(vle.Raoult(P).y(x) - y) for T, x, y in I.VLE_15kPa]
    out(f"- Raoult 모델 vs Table 4: 평균 |ΔT| = {sum(dT)/len(dT):.2f} K, 평균 |Δy| = {sum(dy)/len(dy):.4f}")
    out(f"- 온도: 환류액 {T_cond:.2f} K ({T_cond-273.15:.1f} °C), 탑정 증기 {T_top:.2f} K, 탑저 {T_bot:.2f} K ({T_bot-273.15:.1f} °C)")
    out()

    # 3. Rmin, Nmin, N ----------------------------------------------------
    models = [vle.ConstAlpha(a_exp, f"α 일정 {a_exp:.3f} (Table 4)"), vle.Raoult(P, "Raoult (Table 2)")]
    out("## 3. 최소환류비, 최소단수, 이상단수")
    out("| 평형 모델 | Rmin | R | Nmin Fenske (재비기 제외) | Nmin 작도 (재비기 포함) | 이상단 (재비기 포함) | 공급단 |")
    out("|---|---|---|---|---|---|---|")
    results = {}
    for eq in models:
        Rm, (xp, yp) = mc.R_min(eq, I.xD, I.zF, I.q)
        pinch = mc.check_pinch(eq, Rm * 1.0001, I.xD, I.xB, I.zF, I.q)
        R = I.R_factor * Rm
        a_f = eq.alpha if isinstance(eq, vle.ConstAlpha) else a_raoult_gm
        nf_tot, nf = mc.fenske(a_f, I.xD, I.xB)
        tr = mc.step_stages(eq, I.xD, I.xB, I.zF, I.q, None)
        st = mc.step_stages(eq, I.xD, I.xB, I.zF, I.q, R)
        results[eq.name] = dict(eq=eq, Rm=Rm, R=R, st=st, tr=tr, pinch=pinch, xp=xp, yp=yp)
        out(f"| {eq.name} | {Rm:.3f} | {R:.3f} | {nf:.1f} | {tr['n_stages']:.1f} | {st['n_stages']:.2f} | {st['feed_stage']} |")
    out()
    main_name = models[0].name
    M = results[main_name]
    dN = results[models[1].name]["st"]["n_stages"] - M["st"]["n_stages"]
    out(f"- 보수성 비교: Raoult(Table 2) 결과가 {dN:+.1f}단 (높이 {dN*I.HETP:+.2f} m) → "
        f"**{'Raoult 비교 쪽이 더 보수적' if dN > 0 else '주 계산 쪽이 더 보수적'}**")
    out(f"- 핀치 점검: Rmin에서 (평형선 − 조작선) 최솟값 = {M['pinch']:.2e} → "
        f"{'다른 접점 없음 (q선 교점이 핀치)' if M['pinch'] > -1e-6 else '다른 위치에 접점 있음! Rmin 재검토 필요'}")
    R = M["R"]
    st = M["st"]
    N_total = st["n_stages"]            # 재비기 포함
    N_col = N_total - 1                  # 탑 안의 이론단 (부분 재비기 1단 제외)
    N_rect = st["feed_stage"] - 1        # 정류부 단 수 (공급단은 탈거부에 포함, SP6 p8)
    N_strip = N_col - N_rect
    out(f"- **주 계산: R = {R:.3f}, 이상단 {N_total:.2f} (재비기 포함) → 탑 안 {N_col:.2f}단** "
        f"(정류부 {N_rect}, 탈거부 {N_strip:.2f}, 공급단 = 위에서 {st['feed_stage']}번째)")
    out()

    # 4. 유량과 열량 -------------------------------------------------------
    L = R * D
    V = (R + 1) * D
    Lb = L + I.q * I.F
    Vb = V - (1 - I.q) * I.F
    lamE = lam_EB(T_top)
    lamS = I.DHVAP_ST
    lam_top = I.xD * lamE + (1 - I.xD) * lamS
    lam_bot = I.xB * lamE + (1 - I.xB) * lamS
    Qc = V * lam_top * 1000 / 3600     # kmol/h × kJ/mol → kW
    Qr = Vb * lam_bot * 1000 / 3600
    out("## 4. 내부 유량과 열량")
    out(f"- L = {L:.1f}, V = {V:.1f}, L̄ = {Lb:.1f}, V̄ = {Vb:.1f} kmol/h")
    out(f"- 증발잠열(NIST): EB {lamE:.2f} kJ/mol @ {T_top:.1f} K (열량계 값 보간), 스티렌 {lamS:.1f} kJ/mol (증기압 유도값, 보수적)")
    out(f"- 정몰유량 점검: λ_스티렌/λ_EB = {lamS/lamE:.3f} (차이 {abs(lamS/lamE-1)*100:.1f} %) — 강의 가정(잠열이 거의 같음)의 한계로 기록")
    out(f"- **응축기 열량 Qc = {Qc:.0f} kW,  재비기 열량 Qr = {Qr:.0f} kW**")
    Qr_alt = Vb * (I.xB * lamE + (1 - I.xB) * I.DHVAP_ST_ALT) * 1000 / 3600
    out(f"- (참고) 스티렌 잠열 40.2 kJ/mol 사용 시 Qr = {Qr_alt:.0f} kW")
    out()

    # 5. 범람과 직경 --------------------------------------------------------
    out("## 5. 범람 속도와 탑 직경 (Fig 18.8, Mellapak 250Y)")
    out("| 위치 | T [K] | ρV [kg/m³] | X | Y [m/s] | u_f [m/s] | u 설계 [m/s] | F-factor [Pa^0.5] | 직경 [m] |")
    out("|---|---|---|---|---|---|---|---|---|")
    yB = models[0].y(I.xB)
    sec = {}
    for name, Tv, yv, xl, Lm, Vm in (("탑정(정류부)", T_top, I.xD, I.xD, L, V),
                                     ("탑저(탈거부)", T_bot, yB, I.xB, Lb, Vb)):
        MV = yv * I.M_EB + (1 - yv) * I.M_ST
        ML = xl * I.M_EB + (1 - xl) * I.M_ST
        rhoV = P * MV / (I.R_GAS * Tv)
        rhoL = I.RHO_L
        X = (Lm * ML) / (Vm * MV) * math.sqrt(rhoV / rhoL)
        Y, note = interp_fig18_8(X)
        uf = Y / math.sqrt(rhoV / (rhoL - rhoV))
        u = I.flood_frac * uf
        Qv = Vm * MV / rhoV / 3600          # m3/s
        A = Qv / u
        Dc = math.sqrt(4 * A / math.pi)
        Ff = u * math.sqrt(rhoV)
        sec[name] = dict(Dc=Dc, Ff=Ff, u=u, rhoV=rhoV, X=X, note=note)
        out(f"| {name} | {Tv:.1f} | {rhoV:.3f} | {X:.4f} | {Y:.4f} | {uf:.2f} | {u:.2f} | {Ff:.2f} | {Dc:.2f} |")
    Dc_design = max(s["Dc"] for s in sec.values())
    Dc_round = math.ceil(Dc_design * 10) / 10
    out(f"- X 판독 위치: " + ", ".join(f"{k} {v['note']}" for k, v in sec.items()))
    out(f"- **설계 직경 = {Dc_design:.2f} m → {Dc_round:.1f} m로 올림** (두 위치 중 큰 값)")
    out(f"- ρL = {I.RHO_L} kg/m³ (Wagner 1997, 13.3 kPa 대표값), ρV = PM/RT (이상기체, 예제 21.6 방식)")
    out()

    # 6. 높이 ----------------------------------------------------------------
    Z_rect = N_rect * I.HETP
    Z_strip = N_strip * I.HETP
    Z = Z_rect + Z_strip
    beds_r = math.ceil(Z_rect / I.max_bed_height)
    beds_s = math.ceil(Z_strip / I.max_bed_height)
    out("## 6. 충전 높이 (Z = N × HETP, HETP ≈ H_Oy)")
    out(f"- 정류부 {N_rect}단 × {I.HETP} = {Z_rect:.2f} m,  탈거부 {N_strip:.2f}단 × {I.HETP} = {Z_strip:.2f} m")
    out(f"- **총 충전 높이 Z = {Z:.2f} m**")
    out(f"- 재분배기 기준(≤ {I.max_bed_height:.0f} m/층): 정류부 {beds_r}층, 탈거부 {beds_s}층 → 총 {beds_r+beds_s}층 "
        f"(분배기·지지대·탑저 여유 높이는 강의 근거가 없어 포함하지 않음)")
    out()
    out("| HETP [m] | 충전 높이 Z [m] |")
    out("|---|---|")
    for h in I.HETP_sensitivity:
        mark = " ← 설계값" if abs(h - I.HETP) < 1e-9 else ""
        out(f"| {h} | {N_col*h:.1f}{mark} |")
    out()

    # 7. 압력강하 점검 ---------------------------------------------------------
    dP = I.DP_PER_M * Z / 10.0          # mbar → kPa
    P_bot = P + dP
    T_bot_real = vle.T_bubble(I.xB, P_bot)
    a_bot_real = vle.alpha_raoult(T_bot_real)
    out("## 7. 압력강하 점검 (설계 후)")
    out(f"- ΔP = {I.DP_PER_M} mbar/m × {Z:.2f} m = {dP:.2f} kPa (Salem 실측값, F-factor {I.DP_REF_FFACTOR} 조건)")
    out(f"- 우리 F-factor: " + ", ".join(f"{k} {v['Ff']:.2f}" for k, v in sec.items()) + " Pa^0.5")
    out(f"- 탑저 압력 ≈ {P_bot:.2f} kPa → 탑저 온도 ≈ {T_bot_real:.2f} K (**{T_bot_real-273.15:.1f} °C**), "
        f"기준 {I.T_POLYMER_LIMIT_C:.0f} °C {'이하 ✔' if T_bot_real-273.15 <= I.T_POLYMER_LIMIT_C else '초과 ✖'}")
    out(f"- 한계: 탑저 쪽 실제 α ≈ {a_bot_real:.3f} (Raoult, {P_bot:.1f} kPa) vs 15 kPa 가정 α — "
        f"균일 압력 가정이 탈거부에서 약간 낙관적일 수 있음")
    out()

    # 그림 -------------------------------------------------------------------
    plot_mccabe(M, R, D, "결과_McCabe-Thiele.png", f"McCabe-Thiele (15 kPa, α={a_exp:.3f}, R={R:.2f})")
    plot_mccabe(M, None, D, "결과_전환류.png", "전환류 (최소단수)", total_reflux=True)
    plot_txy(P)
    out("## 그림")
    out("- 결과_McCabe-Thiele.png, 결과_전환류.png, 결과_T-x-y_15kPa.png, 검증_예제21_2.png")

    (OUT / "결과요약.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"\n저장: {OUT / '결과요약.md'}")


def plot_mccabe(M, R, D, fname, title, total_reflux=False):
    eq = M["eq"]
    res = M["tr"] if total_reflux else M["st"]
    fig, ax = plt.subplots(figsize=(7, 7))
    xx = [i / 400 for i in range(401)]
    ax.plot(xx, [eq.y(x) for x in xx], "k-", lw=1.2, label="평형선")
    ax.plot([0, 1], [0, 1], "k:", lw=0.8)
    T, xd, yd = zip(*I.VLE_15kPa)
    ax.plot(xd, yd, "o", ms=3, color="tab:blue", label="Table 4 실측 (15 kPa)")
    if not total_reflux:
        xi, yi = mc.q_line_intersection_with_rect(R, I.xD, I.zF, I.q)
        ax.plot([xi, I.xD], [yi, I.xD], "b-", lw=1, label="정류부 조작선")
        ax.plot([I.xB, xi], [I.xB, yi], "g-", lw=1, label="탈거부 조작선")
        ax.plot([I.zF, I.zF], [I.zF, yi], "m--", lw=1, label="q선")
    px, py = zip(*res["points"])
    ax.plot(px, py, "r-", lw=0.6, label=f"계단 {res['n_stages']:.1f}단")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("x (EB 몰분율)"); ax.set_ylabel("y (EB 몰분율)")
    ax.set_title(title); ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(OUT / fname, dpi=150); plt.close(fig)


def plot_txy(P):
    fig, ax = plt.subplots(figsize=(7, 5))
    xx = [i / 200 for i in range(201)]
    Tb = [vle.T_bubble(x, P) for x in xx]
    Td = [vle.T_dew(x, P) for x in xx]
    ax.plot(xx, Tb, "b-", label="기포점 (Raoult, Table 2)")
    ax.plot(xx, Td, "r-", label="이슬점 (Raoult, Table 2)")
    T, xd, yd = zip(*I.VLE_15kPa)
    ax.plot(xd, T, "bo", ms=4, label="x 실측 (Table 4)")
    ax.plot(yd, T, "rs", ms=4, label="y 실측 (Table 4)")
    ax.set_xlabel("x, y (EB 몰분율)"); ax.set_ylabel("T [K]")
    ax.set_title(f"T-x-y 선도 ({P} kPa)"); ax.legend()
    fig.tight_layout(); fig.savefig(OUT / "결과_T-x-y_15kPa.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    main()
