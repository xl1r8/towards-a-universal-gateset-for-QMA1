from sage.all import *

K = QQbar

idm = lambda n: identity_matrix(K, n)
I2 = idm(2)
I4 = idm(4)
dagger = lambda x: x.conjugate_transpose()
td = 12
qd = 4
qbra = lambda i: matrix(K, [[0] * i + [1] + [0] * (2 - i - 1)])
qket = lambda i: dagger(qbra(i))
tbra = lambda i: matrix(K, [[0] * i + [1] + [0] * (td - i - 1)])
tket = lambda i: dagger(tbra(i))
qkb = lambda a, b: qket(a) * qbra(b)
tkb = lambda a, b: tket(a) * tbra(b)


def kb(ketx):
    return ketx * dagger(ketx)


U_trans = lambda t0, t1, U: tp(tkb(t0, t0) + tkb(t1, t1), idm(qd)) - tp(tkb(t0, t1), dagger(U)) - tp(tkb(t1, t0), U)
I_trans = lambda t0, t1: U_trans(t0, t1, idm(qd))


def tp(*args):
    return reduce(lambda a, b: a.tensor_product(b), args)


def split_gadget(ti, t0, t1):  # ti: input, t0: 0 in control, t1: 1 in control
    T = tkb(ti, ti) - tkb(ti, t0) - tkb(ti, t1) + \
        tkb(t0, t0) - tkb(t0, ti) + tkb(t0, t1) + \
        tkb(t1, t1) - tkb(t1, ti) + tkb(t1, t0)
    return T


def setup_dim(t, q):
    global td, qd
    td = t
    qd = q


def test_split_gadget():
    setup_dim(3, 2)
    pen = lambda t, b: tp(tkb(t, t), kb(qket(b)))
    Hpen = pen(1, 1) + pen(2, 0)
    Hsplit = tp(split_gadget(0, 2, 1), I2)
    H = Hpen + Hsplit

    assert H.nrows() - H.rank() == 2
    assert H.is_hermitian()
    assert H.is_positive_semidefinite()

    tq = lambda t, b: tp(tket(t), qket(b))
    v0 = tq(0, 0) + tq(1, 0)
    v1 = tq(0, 1) + tq(2, 1)
    assert H * v0 == 0
    assert H * v1 == 0


def test_X_gadget():
    setup_dim(6, 2)

    B = [qket(0) + qket(1), qket(0) - qket(1)]
    pen = lambda t, b: tp(tkb(t, t), kb(B[b]))
    Hpen = (pen(1, 1) +
            pen(2, 0) +

            pen(3, 1) +
            pen(4, 0))

    Hsplit = tp(split_gadget(0, 1, 2) + split_gadget(5, 3, 4), I2)

    Htrans = U_trans(2, 4, -I2) + I_trans(1, 3)

    H = Hpen + Hsplit + Htrans

    tq = lambda t, b: tp(tket(t), B[b])
    v0 = tq(0, 0) + tq(1, 0) + tq(3, 0) + tq(5, 0)
    v1 = tq(0, 1) + tq(2, 1) - tq(4, 1) - tq(5, 1)
    assert H * v0 == 0
    assert H * v1 == 0

    assert H.is_hermitian()
    assert H.is_positive_semidefinite()
    assert H.nrows() - H.rank() == 2


def test_T_gadget():
    setup_dim(6, 2)

    w = exp(I * pi / 4)

    pen = lambda t, b: tp(tkb(t, t), kb(qket(b)))
    Hpen = (pen(1, 1) +
            pen(2, 0) +

            pen(3, 1) +
            pen(4, 0))

    Hsplit = tp(split_gadget(0, 1, 2) + split_gadget(5, 3, 4), I2)

    Htrans = U_trans(2, 4, w * I2) + I_trans(1, 3)

    H = Hpen + Hsplit + Htrans

    tq = lambda t, b: tp(tket(t), qket(b))
    v0 = tq(0, 0) + tq(1, 0) + tq(3, 0) + tq(5, 0)
    v1 = tq(0, 1) + tq(2, 1) + w * tq(4, 1) + w * tq(5, 1)
    assert H * v0 == 0
    assert H * v1 == 0

    assert H.is_hermitian()
    assert H.is_positive_semidefinite()
    assert H.nrows() - H.rank() == 2


def test_H_gadget():
    setup_dim(6, 2)

    H_gate = matrix(K, [[1, 1], [1, -1]]) / sqrt(2)

    a = K(sqrt(2 + sqrt(2)) / 2)
    b = K(sqrt(2 - sqrt(2)) / 2)
    # eigenstates of H
    B = [a * qket(0) + b * qket(1), -b * qket(0) + a * qket(1)]

    assert H_gate * B[0] == B[0]
    assert H_gate * B[1] == -B[1]

    pen = lambda t, b: tp(tkb(t, t), kb(B[b]))
    Hpen = (pen(1, 1) +
            pen(2, 0) +

            pen(3, 1) +
            pen(4, 0))

    Hsplit = tp(split_gadget(0, 1, 2) + split_gadget(5, 3, 4), I2)

    Htrans = U_trans(2, 4, -I2) + I_trans(1, 3)

    H = Hpen + Hsplit + Htrans

    tq = lambda t, b: tp(tket(t), B[b])
    v0 = tq(0, 0) + tq(1, 0) + tq(3, 0) + tq(5, 0)
    v1 = tq(0, 1) + tq(2, 1) - tq(4, 1) - tq(5, 1)
    assert H * v0 == 0
    assert H * v1 == 0

    assert H.is_hermitian()
    assert H.is_positive_semidefinite()
    assert H.nrows() - H.rank() == 2

    k = CyclotomicField(8)
    matrix(k, H)  # check that H is in QQ(zeta_8)


def test_CNOT_gadget():
    setup_dim(12, 4)

    B = [qket(0) + qket(1), qket(0) - qket(1)]
    pen = lambda t, i, b: tp(tkb(t, t),
                             tp(qkb(b, b), I2) if i == 0 else tp(I2, kb(B[b])))
    Hpen = (pen(1, 0, 1) +
            pen(2, 0, 0) +

            pen(3, 0, 1) +
            pen(4, 0, 0) + pen(4, 1, 1) +
            pen(5, 0, 0) + pen(5, 1, 0) +

            pen(6, 0, 1) +
            pen(7, 0, 0) + pen(7, 1, 1) +
            pen(8, 0, 0) + pen(8, 1, 0) +

            pen(9, 0, 1) +
            pen(10, 0, 0))

    Hsplit = tp(split_gadget(0, 1, 2) + split_gadget(2, 4, 5) +
                split_gadget(10, 7, 8) + split_gadget(11, 9, 10), I4)

    Htrans = (U_trans(5, 8, -I4) + I_trans(4, 7) +
              I_trans(3, 6) + I_trans(1, 3) + I_trans(6, 9))

    H = Hpen + Hsplit + Htrans

    tq = lambda t, b0, b1: tp(tket(t), qket(b0), B[b1])
    v00 = tq(0, 0, 0) + tq(1, 0, 0) + tq(3, 0, 0) + tq(6, 0, 0) + tq(9, 0, 0) + tq(11, 0, 0)
    v01 = tq(0, 0, 1) + tq(1, 0, 1) + tq(3, 0, 1) + tq(6, 0, 1) + tq(9, 0, 1) + tq(11, 0, 1)
    v10 = tq(0, 1, 0) + tq(2, 1, 0) + tq(4, 1, 0) + tq(7, 1, 0) + tq(10, 1, 0) + tq(11, 1, 0)
    v11 = tq(0, 1, 1) + tq(2, 1, 1) + tq(5, 1, 1) - tq(8, 1, 1) - tq(10, 1, 1) - tq(11, 1, 1)

    # verify that H is Hermitian, PSD, and its nullspace is Span{v00, v01, v10, v11}
    assert H * v00 == 0
    assert H * v01 == 0
    assert H * v10 == 0
    assert H * v11 == 0

    assert H.is_hermitian()
    assert H.is_positive_semidefinite()
    assert H.nrows() - H.rank() == 4


def test_controlled_S_gadget():
    setup_dim(12, 4)

    pen = lambda t, i, b: tp(tkb(t, t),
                             tp(qkb(b, b), I2) if i == 0 else tp(I2, qkb(b, b)))
    Hpen = (pen(1, 0, 1) +
            pen(2, 0, 0) +

            pen(3, 0, 1) +
            pen(4, 0, 0) + pen(4, 1, 1) +
            pen(5, 0, 0) + pen(5, 1, 0) +

            pen(6, 0, 1) +
            pen(7, 0, 0) + pen(7, 1, 1) +
            pen(8, 0, 0) + pen(8, 1, 0) +

            pen(9, 0, 1) +
            pen(10, 0, 0))

    Hsplit = tp(split_gadget(0, 1, 2) + split_gadget(2, 4, 5) +
                split_gadget(10, 7, 8) + split_gadget(11, 9, 10), I4)

    Htrans = (U_trans(5, 8, I4 * I) + I_trans(4, 7) +
              I_trans(3, 6) + I_trans(1, 3) + I_trans(6, 9))

    H = Hpen + Hsplit + Htrans

    tq = lambda t, b0, b1: tp(tket(t), qket(b0), qket(b1))
    v00 = tq(0, 0, 0) + tq(1, 0, 0) + tq(3, 0, 0) + tq(6, 0, 0) + tq(9, 0, 0) + tq(11, 0, 0)
    v01 = tq(0, 0, 1) + tq(1, 0, 1) + tq(3, 0, 1) + tq(6, 0, 1) + tq(9, 0, 1) + tq(11, 0, 1)
    v10 = tq(0, 1, 0) + tq(2, 1, 0) + tq(4, 1, 0) + tq(7, 1, 0) + tq(10, 1, 0) + tq(11, 1, 0)
    v11 = tq(0, 1, 1) + tq(2, 1, 1) + tq(5, 1, 1) + I * tq(8, 1, 1) + I * tq(10, 1, 1) + I * tq(11, 1, 1)

    # verify that H is Hermitian, PSD, and its nullspace is Span{v00, v01, v10, v11}
    assert H * v00 == 0
    assert H * v01 == 0
    assert H * v10 == 0
    assert H * v11 == 0

    assert H.is_hermitian()
    assert H.is_positive_semidefinite()
    assert H.nrows() - H.rank() == 4


def test_H_eigenbasis():
    # see also https://meirizarrygelpi.github.io/posts/physics/hadamard-eigen-basis/index.html
    H_gate = matrix(K, [[1, 1], [1, -1]]) / sqrt(2)
    a = sqrt(2 + sqrt(2)) / 2
    b = sqrt(2 - sqrt(2)) / 2
    H_plus = a * qket(0) + b * qket(1)
    H_minus = -b * qket(0) + a * qket(1)
    assert H_gate == kb(H_plus) - kb(H_minus)
    assert H_gate * H_plus == H_plus
    assert H_gate * H_minus == -H_minus
    assert a == cos(pi / 8)
    assert b == sin(pi / 8)
    k = CyclotomicField(8)
    matrix(k, kb(H_plus))
    matrix(k, kb(H_minus))
    assert a * a == sqrt(2) / 4 + 1 / 2
    assert b * b == -sqrt(2) / 4 + 1 / 2
    assert a * b == sqrt(2) / 4


test_split_gadget()
test_X_gadget()
test_T_gadget()
test_H_gadget()
test_CNOT_gadget()
test_controlled_S_gadget()
test_H_eigenbasis()
