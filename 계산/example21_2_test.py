"""
방법 검증 (Method Implementation and Test)
강의 SP6 p30~p36, Example 21.2 (벤젠–톨루엔, α ≈ 2.5, R = 3.5, 포화액 공급)를
같은 함수로 다시 풀어 교재 답과 비교한다.

교재 답 (SP6 p32, p33 그림):
  D = 153.4 kmol/h, B = 196.6 kmol/h
  이상단 11개 (재비기 포함), 원료 공급단 7번째
※ 교재는 실측 평형곡선(Fig. 21.14)을 썼고 여기서는 α = 2.5 일정이므로 약간 차이가 날 수 있다.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from vle import ConstAlpha
import mccabe as mc

OUT = Path(__file__).parent / "결과"
OUT.mkdir(exist_ok=True)
plt.rcParams["font.family"] = "Malgun Gothic"
plt.rcParams["axes.unicode_minus"] = False


def run(verbose=True):
    F, xF, xD, xB, q, R = 350.0, 0.440, 0.974, 0.0235, 1.0, 3.5
    eq = ConstAlpha(2.5)
    D, B = mc.mass_balance(F, xF, xD, xB)
    Rm, _ = mc.R_min(eq, xD, xF, q)
    res = mc.step_stages(eq, xD, xB, xF, q, R)
    ok = (abs(D - 153.4) < 0.5) and (abs(res["n_int"] - 11) <= 1) and (abs(res["feed_stage"] - 7) <= 1)
    if verbose:
        print("=" * 64)
        print("[방법 검증] Example 21.2 (벤젠-톨루엔)")
        print(f"  D = {D:.1f} kmol/h (교재 153.4),  B = {B:.1f} kmol/h (교재 196.6)")
        print(f"  Rmin (α=2.5 일정) = {Rm:.3f}")
        print(f"  이상단 수 (재비기 포함) = {res['n_stages']:.2f}  → 계단 {res['n_int']}개 (교재 11)")
        print(f"  원료 공급단 = {res['feed_stage']}번째 (교재 7)")
        print(f"  판정: {'±1단 이내 일치' if ok else '확인 필요'} (교재는 실측 평형곡선, 여기는 α=2.5 일정)")
    # 그림
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    xx = [i / 200 for i in range(201)]
    ax.plot(xx, [eq.y(x) for x in xx], "k-", label="평형선 (α=2.5)")
    ax.plot([0, 1], [0, 1], "k:", lw=0.8)
    px, py = zip(*res["points"])
    ax.plot(px, py, "r-", lw=1, label="계단 작도")
    ax.axvline(xF, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("x (벤젠)"); ax.set_ylabel("y (벤젠)")
    ax.set_title(f"Example 21.2 재현: {res['n_int']}단, 공급단 {res['feed_stage']}")
    ax.legend(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.tight_layout(); fig.savefig(OUT / "검증_예제21_2.png", dpi=150); plt.close(fig)
    return dict(D=D, B=B, Rmin=Rm, n=res["n_stages"], n_int=res["n_int"], feed=res["feed_stage"], ok=ok)


if __name__ == "__main__":
    run()
