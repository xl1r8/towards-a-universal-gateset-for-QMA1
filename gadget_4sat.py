from sage.all import *

K = QQbar

idm = lambda n: identity_matrix(K, n)
I2 = idm(2)
I4 = idm(4)
dagger = lambda x: x.conjugate_transpose()
td = 3
qbra = lambda i: matrix(K, [[0] * i + [1] + [0] * (2 - i - 1)])
qket = lambda i: dagger(qbra(i))
tbra = lambda i: matrix(K, [[0] * i + [1] + [0] * (td - i - 1)])
tket = lambda i: dagger(tbra(i))
qkb = lambda a, b: qket(a) * qbra(b)
tkb = lambda a, b: tket(a) * tbra(b)


def kb(ketx):
    return ketx * dagger(ketx)


def I_trans(t0, t1):
    return kb(tket(t0) - tket(t1))


def tp(*args):
    return reduce(lambda a, b: a.tensor_product(b), args)


def test_X_gadget():
    tq = lambda t, b: tp(tket(t), qket(b))
    X_trans = kb(tq(0, 0) - tq(1, 1)) + kb(tq(0, 1) - tq(1, 0))
    H = X_trans + tp(I_trans(1, 2), I2)

    v0 = tq(0, 0) + tq(1, 1) + tq(2, 1)
    v1 = tq(0, 1) + tq(1, 0) + tq(2, 0)
    assert H * v0 == 0
    assert H * v1 == 0

    assert H.nrows() - H.rank() == 2
    assert H.is_hermitian()
    assert H.is_positive_semidefinite()


def test_CX_gadget():
    tq = lambda t, b0, b1: tp(tket(t), qket(b0), qket(b1))
    X_trans = (
            kb(tq(0, 0, 0) - tq(1, 0, 0)) +
            kb(tq(0, 0, 1) - tq(1, 0, 1)) +
            kb(tq(0, 1, 0) - tq(1, 1, 1)) +
            kb(tq(0, 1, 1) - tq(1, 1, 0))
    )
    H = X_trans + tp(I_trans(1, 2), I4)

    v00 = tq(0, 0, 0) + tq(1, 0, 0) + tq(2, 0, 0)
    v01 = tq(0, 0, 1) + tq(1, 0, 1) + tq(2, 0, 1)
    v10 = tq(0, 1, 0) + tq(1, 1, 1) + tq(2, 1, 1)
    v11 = tq(0, 1, 1) + tq(1, 1, 0) + tq(2, 1, 0)
    assert H * v00 == 0
    assert H * v01 == 0
    assert H * v10 == 0
    assert H * v11 == 0

    assert H.nrows() - H.rank() == 4
    assert H.is_hermitian()
    assert H.is_positive_semidefinite()


def test_CCX_gadget():
    tq = lambda t, b0, b1: tp(tket(t), I2, qket(b0), qket(b1))
    CX_trans = (
            kb(tq(0, 0, 0) - tq(1, 0, 0)) +
            kb(tq(0, 0, 1) - tq(1, 0, 1)) +
            kb(tq(0, 1, 0) - tq(1, 1, 1)) +
            kb(tq(0, 1, 1) - tq(1, 1, 0))
    )
    H = CX_trans + tp(I_trans(1, 2), qkb(1, 1), I2, I2) + tp(I_trans(0, 2), qkb(0, 0), I2, I2)

    tq = lambda t, b0, b1, b2: tp(tket(t), qket(b0), qket(b1), qket(b2))
    v000 = tq(0, 0, 0, 0) + tq(1, 0, 0, 0) + tq(2, 0, 0, 0)
    v001 = tq(0, 0, 0, 1) + tq(1, 0, 0, 1) + tq(2, 0, 0, 1)
    v010 = tq(0, 0, 1, 0) + tq(1, 0, 1, 1) + tq(2, 0, 1, 0)
    v011 = tq(0, 0, 1, 1) + tq(1, 0, 1, 0) + tq(2, 0, 1, 1)
    v100 = tq(0, 1, 0, 0) + tq(1, 1, 0, 0) + tq(2, 1, 0, 0)
    v101 = tq(0, 1, 0, 1) + tq(1, 1, 0, 1) + tq(2, 1, 0, 1)
    v110 = tq(0, 1, 1, 0) + tq(1, 1, 1, 1) + tq(2, 1, 1, 1)
    v111 = tq(0, 1, 1, 1) + tq(1, 1, 1, 0) + tq(2, 1, 1, 0)

    assert H * v000 == 0
    assert H * v001 == 0
    assert H * v010 == 0
    assert H * v011 == 0
    assert H * v100 == 0
    assert H * v101 == 0
    assert H * v110 == 0
    assert H * v111 == 0

    assert H.nrows() - H.rank() == 8
    assert H.is_hermitian()
    assert H.is_positive_semidefinite()


def test_HH_gadget():
    B = [qket(0) + qket(1), qket(0) - qket(1)]
    tq = lambda t, b: tp(tket(t), qket(b), I2)
    tqH = lambda t, b: tp(tket(t), B[b], I2)
    # sqrt(2)*H
    # |00> - |10> - |11>
    # |01> - |10> + |11>
    H_trans = (
            kb(tq(0, 0) - tqH(1, 0)) +
            kb(tq(0, 1) - tqH(1, 1))
    )
    tq = lambda t, b: tp(tket(t), I2, qket(b))
    tqH = lambda t, b: tp(tket(t), I2, B[b])
    # sqrt(1/2)*H
    # |10> + |11> - |00> = -(|00> - |10> - |11>)
    # |10> - |11> - |01> = -(|01> - |10> + |11>)
    H_trans += (
            kb(tqH(1, 0) - tq(2, 0)) +
            kb(tqH(1, 1) - tq(2, 1))
    )
    # alternative:
    # H_trans += (
    #         kb(2*tq(1, 0) - tqH(2, 0)) +
    #         kb(2*tq(1, 1) - tqH(2, 1))
    # )
    H = H_trans

    v00 = tp(tket(0), qket(0), qket(0)) + tp(tket(1), B[0], qket(0)) / 2 + tp(tket(2), B[0], B[0]) / 2
    v01 = tp(tket(0), qket(0), qket(1)) + tp(tket(1), B[0], qket(1)) / 2 + tp(tket(2), B[0], B[1]) / 2
    v10 = tp(tket(0), qket(1), qket(0)) + tp(tket(1), B[1], qket(0)) / 2 + tp(tket(2), B[1], B[0]) / 2
    v11 = tp(tket(0), qket(1), qket(1)) + tp(tket(1), B[1], qket(1)) / 2 + tp(tket(2), B[1], B[1]) / 2
    assert H * v00 == 0
    assert H * v01 == 0
    assert H * v10 == 0
    assert H * v11 == 0

    assert H.nrows() - H.rank() == 4
    assert H.is_hermitian()
    assert H.is_positive_semidefinite()


test_X_gadget()
test_CX_gadget()
test_CCX_gadget()
test_HH_gadget()
