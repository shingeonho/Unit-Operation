"""
기액평형(VLE) 모델 — 강의 SP6의 상대휘발도 개념만 사용
  - ConstAlpha : y = αx / (1 + (α-1)x)            (SP6 p4, p25)
  - Raoult     : α = P_EB'/P_ST' (Raoult 법칙),  P'는 논문 Table 2 실측값 보간
"""
import math
import inputs as I


def _interp_lnP(table, T):
    """Table 2 (T, P) 사이를 ln P – 1/T 직선으로 보간 (범위 밖은 끝 구간 연장)."""
    xs = [1.0 / t for t, _ in table]       # T 증가 → 1/T 감소
    ys = [math.log(p) for _, p in table]
    x = 1.0 / T
    # xs는 감소 순서. 해당 구간 찾기
    if x >= xs[0]:
        i = 0
    elif x <= xs[-1]:
        i = len(xs) - 2
    else:
        i = next(k for k in range(len(xs) - 1) if xs[k] >= x >= xs[k + 1])
    x0, x1, y0, y1 = xs[i], xs[i + 1], ys[i], ys[i + 1]
    return math.exp(y0 + (y1 - y0) * (x - x0) / (x1 - x0))


def psat_EB(T):
    return _interp_lnP(I.PSAT_EB, T)


def psat_ST(T):
    return _interp_lnP(I.PSAT_ST, T)


def _bisect(f, a, b, tol=1e-9, it=200):
    fa = f(a)
    for _ in range(it):
        m = 0.5 * (a + b)
        fm = f(m)
        if fa * fm <= 0:
            b = m
        else:
            a, fa = m, fm
        if b - a < tol:
            break
    return 0.5 * (a + b)


def T_bubble(x, P):
    """기포점 [K]: x·P_EB' + (1-x)·P_ST' = P"""
    return _bisect(lambda T: x * psat_EB(T) + (1 - x) * psat_ST(T) - P, 280.0, 450.0)


def T_dew(y, P):
    """이슬점 [K]: y·P/P_EB' + (1-y)·P/P_ST' = 1"""
    return _bisect(lambda T: 1.0 - (y * P / psat_EB(T) + (1 - y) * P / psat_ST(T)), 280.0, 450.0)


def alpha_raoult(T):
    return psat_EB(T) / psat_ST(T)


class ConstAlpha:
    """상대휘발도 일정 모델 (강의 SP6 p25)"""
    def __init__(self, alpha, name="const-alpha"):
        self.alpha = alpha
        self.name = name

    def y(self, x):
        a = self.alpha
        return a * x / (1 + (a - 1) * x)

    def x(self, y):
        a = self.alpha
        return y / (a - (a - 1) * y)


class Raoult:
    """압력 P에서 Raoult 법칙 + Table 2 증기압 (온도에 따라 α가 변함)"""
    def __init__(self, P, name="Raoult(Table 2)"):
        self.P = P
        self.name = name

    def y(self, x):
        T = T_bubble(x, self.P)
        return x * psat_EB(T) / self.P

    def x(self, y):
        T = T_dew(y, self.P)
        return y * self.P / psat_EB(T)


def alpha_from_data(data):
    """실측 (T, x, y) 각 점의 α와 기하평균"""
    al = [(y / (1 - y)) / (x / (1 - x)) for _, x, y in data]
    gm = math.exp(sum(math.log(a) for a in al) / len(al))
    return al, gm
