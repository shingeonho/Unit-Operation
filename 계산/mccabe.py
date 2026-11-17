"""
McCabe-Thiele 계산 함수 — 강의 SP6 p9~p28의 식만 사용
  물질수지, 조작선, q선, 최소환류비, Fenske, 단 작도(stepping)
eq 객체는 eq.y(x), eq.x(y)를 가진 평형 모델 (vle.py)
"""
import math


def mass_balance(F, zF, xD, xB):
    """SP6 p9: D/F = (zF-xB)/(xD-xB)"""
    D = F * (zF - xB) / (xD - xB)
    return D, F - D


def q_line_intersection_with_rect(R, xD, zF, q):
    """정류부 조작선과 q선의 교점 (SP6 p21)"""
    if abs(q - 1.0) < 1e-12:
        x = zF
    else:
        # y = R/(R+1) x + xD/(R+1)  와  y = q/(q-1) x - zF/(q-1)
        a1, b1 = R / (R + 1), xD / (R + 1)
        a2, b2 = q / (q - 1), -zF / (q - 1)
        x = (b2 - b1) / (a1 - a2)
    y = R / (R + 1) * x + xD / (R + 1)
    return x, y


def q_line_intersection_with_eq(eq, zF, q):
    """q선과 평형선의 교점 (x', y')"""
    if abs(q - 1.0) < 1e-12:
        return zF, eq.y(zF)
    f = lambda x: eq.y(x) - (q / (q - 1) * x - zF / (q - 1))
    a, b = 1e-6, 1 - 1e-6
    fa = f(a)
    for _ in range(200):
        m = 0.5 * (a + b)
        if fa * f(m) <= 0:
            b = m
        else:
            a, fa = m, f(m)
    x = 0.5 * (a + b)
    return x, eq.y(x)


def R_min(eq, xD, zF, q):
    """SP6 p27: R_Dm = (xD - y')/(y' - x')"""
    xp, yp = q_line_intersection_with_eq(eq, zF, q)
    return (xD - yp) / (yp - xp), (xp, yp)


def check_pinch(eq, R, xD, xB, zF, q, n=400):
    """최소환류비에서 조작선이 평형선을 다른 곳에서 넘지 않는지 확인 (SP6 p28)
    반환: (평형선-조작선)의 최솟값 (음수면 접선 pinch가 다른 곳에 있음)"""
    xi, yi = q_line_intersection_with_rect(R, xD, zF, q)
    worst = 1.0
    for k in range(n + 1):
        x = xB + (xD - xB) * k / n
        if x >= xi:
            yop = R / (R + 1) * x + xD / (R + 1)
        else:
            yop = xB + (yi - xB) / (xi - xB) * (x - xB)
        worst = min(worst, eq.y(x) - yop)
    return worst


def fenske(alpha, xD, xB):
    """SP6 p26: (총 이상단, 재비기 제외한 단수 N_min)"""
    n_total = math.log(xD * (1 - xB) / (xB * (1 - xD))) / math.log(alpha)
    return n_total, n_total - 1


def step_stages(eq, xD, xB, zF, q, R=None, max_steps=1000):
    """McCabe-Thiele 단 작도.  R=None 이면 전환류(조작선 y=x).
    전응축기: (xD, xD)에서 시작.  q선 교점의 x보다 작아지면 탈거부 조작선으로 전환(SP6 p23).
    반환 dict: points(계단 꼭짓점), n_stages(재비기 포함, 소수), feed_stage, xs(각 단 액체 조성)"""
    if R is None:
        op_rect = op_strip = lambda x: x
        xi = None
    else:
        xi, yi = q_line_intersection_with_rect(R, xD, zF, q)
        op_rect = lambda x: R / (R + 1) * x + xD / (R + 1)
        op_strip = lambda x: xB + (yi - xB) / (xi - xB) * (x - xB)

    pts = [(xD, xD)]
    xs = []
    y = xD
    stripping = False
    feed_stage = None
    x_prev = xD
    for n in range(1, max_steps + 1):
        x = eq.x(y)
        pts.append((x, y))
        if x <= xB:
            frac = (x_prev - xB) / (x_prev - x)
            xs.append(x)
            return dict(points=pts, n_stages=n - 1 + frac, n_int=n, feed_stage=feed_stage, xs=xs)
        xs.append(x)
        if xi is not None and not stripping and x < xi:
            stripping = True
            feed_stage = n
        y_next = (op_strip if stripping else op_rect)(x)
        pts.append((x, y_next))
        x_prev, y = x, y_next
    raise RuntimeError("단 수가 max_steps를 넘음 (R이 Rmin에 너무 가까움?)")
