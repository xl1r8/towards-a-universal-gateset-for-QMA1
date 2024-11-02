from sage.all import *
import itertools

v0 = 'v0'
vertex_names = 'xx, a2, a3, a4, b2, b3, b4'.split(', ')


def make_graph(n):
    """ Constructs graph G_n """
    names = x, a2, a3, a4, b2, b3, b4 = [[f'{i}.{v}' for i in range(n)] for v in vertex_names]
    qubit_vertices = [[vv[i] for vv in names] for i in range(n)]
    G = Graph()
    for i in range(n):
        G.add_edges([[x[i], a3[i]], [a3[i], a2[i]], [a2[i], a4[i]], [a4[i], x[i]], [x[i], b3[i]], [b3[i], b2[i]],
                     [b2[i], b4[i]], [b4[i], x[i]]])
        for j in range(i):
            G.add_edges(itertools.product(qubit_vertices[i], qubit_vertices[j]))
    return G, names


def thicken(G):
    # apply Lemma 7.14 to G_cycle
    G_thicken = Graph()
    for (u, v) in G.edges(labels=False, sort_vertices=True):
        G_thicken.add_edges([(f'{u}.0', f'{v}.0'),
                             (f'{u}.1', f'{v}.1'),
                             (f'{u}.0', f'{v}.1'), ])
    return G_thicken


def fill_cycle(G_base, G_cycle, f={}):  # assume cycle is subgraph of base
    merge = {v: [] for v in G_base}
    G_thicken = thicken(G_cycle)
    for v in G_cycle:
        G_thicken.add_edges([(f'{v}.{0}', f'{v}.{1}'),
                             (f'{v}.{1}', v0)])
        merge[f.get(v, v)].append(f'{v}.{0}')
    # G_thicken.plot().save("G_thicken.png")
    G_fill = G_thicken.union(G_base)
    # G_thicken.plot().save("G_fill1.png")
    for v, vv in merge.items():
        G_fill.merge_vertices([v] + vv)
    # G_fill.plot().save("G_fill2.png")
    return G_fill


def cell2vec(cell2dim, cells):
    v = zero_vector(QQ, len(cell2dim))
    for cell in cells:
        c, s = sort_and_sign(cell)
        d, s2 = cell2dim[c]
        v[d] += s * s2
    return v


def sort_and_sign(a):
    s = tuple(sorted(a))
    perm = Permutation([s.index(x) + 1 for x in a])
    return s, perm.sign()


def demo_fill_cycle():
    G, names = make_graph(1)
    G.plot().save("G.png")
    x, a2, a3, a4, b2, b3, b4 = names
    G_cycle = Graph()
    G_cycle.add_edges([(x[0], a3[0]), (a3[0], a2[0]), (a2[0], a4[0]), (a4[0], x[0]), ])
    print(G_cycle.clique_complex().homology())
    G_cycle.plot().save("G_cycle.png")
    G_fill = fill_cycle(G, G_cycle)
    G_fill.plot(layout="spring").save("G_fill.png")

    S = G_fill.clique_complex()
    print(S.homology())
    C = S.chain_complex(base_ring=QQ)
    d = C.differential()
    print(d)
    L = d[1].transpose() * d[1] + d[2] * d[2].transpose()
    print(L.right_kernel())

    cell2dim = {}
    for i, cell in enumerate(S._n_cells_sorted(1)):
        c, s = sort_and_sign(cell)
        cell2dim[c] = (i, s)
    print(cell2dim)
    ket0 = cell2vec(cell2dim, [(x[0], a3[0]), (a3[0], a2[0]), (a2[0], a4[0]), (a4[0], x[0]), ])
    ket1 = cell2vec(cell2dim, [(x[0], b3[0]), (b3[0], b2[0]), (b2[0], b4[0]), (b4[0], x[0]), ])
    print(ket1)
    print(L * ket1 == 0)

    V = d[1].dense_matrix().right_kernel()
    W = d[2].dense_matrix().transpose().image()
    print(f'dim(V) = {V.dimension()}, dim(W) = {W.dimension()}')
    Q = V.quotient(W)
    print(f'{ket0 in V = }, {ket0 in W = }')
    print(f'{ket1 in V = }, {ket1 in W = }')


def state_0m1():
    G, names = make_graph(1)
    x, a2, a3, a4, b2, b3, b4 = names
    x1 = '0.x1'
    G_cycle = Graph()
    # G_cycle.add_edges(  # as in Figure 7
    #     [(a2[0], a3[0]), (a3[0], x[0]), (x[0], b4[0]), (b4[0], b2[0]), (b2[0], b3[0]), (b3[0], x1), (x1, a4[0]),
    #      (a4[0], a2[0]), ])
    G_cycle.add_edges(
        [(a2[0], a3[0]), (a3[0], x[0]), (x[0], b3[0]), (b3[0], b2[0]), (b2[0], b4[0]), (b4[0], x1),
         (x1, a4[0]), (a4[0], a2[0]), ])
    G_fill = fill_cycle(G, G_cycle, {x1: x[0]})

    S = G_fill.clique_complex()
    assert str(S.homology()) == '{0: 0, 1: Z, 2: 0}'

    C = S.chain_complex(base_ring=QQ)
    d = C.differential()

    cell2dim = {}
    for i, cell in enumerate(S._n_cells_sorted(1)):
        c, s = sort_and_sign(cell)
        cell2dim[c] = (i, s)

    ket1 = cell2vec(cell2dim, [(x[0], b3[0]), (b3[0], b2[0]), (b2[0], b4[0]), (b4[0], x[0]), ])
    ket0 = cell2vec(cell2dim, [(x[0], a3[0]), (a3[0], a2[0]), (a2[0], a4[0]), (a4[0], x[0]), ])
    ket_plus = ket0 + ket1
    ket_minus = ket0 - ket1

    V = d[1].dense_matrix().right_kernel()
    W = d[2].dense_matrix().transpose().image()
    Q = V.quotient(W)
    assert Q.dimension() == 1

    assert ket_minus in V and ket_minus in W
    assert ket_plus in V and ket_plus not in W


def demo_fill_cycle_0_1_1():
    G, names = make_graph(1)
    G.plot().save("G.png")
    x, a2, a3, a4, b2, b3, b4 = names
    x1 = 'x1'
    x2 = 'x2'
    b2_ = '0.b2_'
    b3_ = '0.b3_'
    b4_ = '0.b4_'
    G_cycle = Graph()
    G_cycle.add_edges([
        (a2[0], a3[0]),
        (a3[0], x[0]),
        (x[0], b3[0]),
        (b3[0], b2[0]),
        (b2[0], b4[0]),
        (b4[0], x2),
        (x2, b3_),
        (b3_, b2_),
        (b2_, b4_),
        (b4_, x1),
        (x1, a4[0]),
        (a4[0], a2[0]),
    ])
    G_cycle.plot().save("G_cycle.png")
    G_fill = fill_cycle(G, G_cycle, {x1: x[0], x2: x[0], b2_: b2[0], b3_: b3[0], b4_: b4[0]})

    S = G_fill.clique_complex()
    print(S.homology())
    C = S.chain_complex(base_ring=QQ)
    d = C.differential()
    print(d)
    L = d[1].transpose() * d[1] + d[2] * d[2].transpose()
    L_ker = L.right_kernel()
    print(L_ker)

    cell2dim = {}
    for i, cell in enumerate(S._n_cells_sorted(1)):
        c, s = sort_and_sign(cell)
        cell2dim[c] = (i, s)
    print(cell2dim)
    ket1 = cell2vec(cell2dim, [(x[0], b3[0]), (b3[0], b2[0]), (b2[0], b4[0]), (b4[0], x[0]), ])
    ket0 = cell2vec(cell2dim, [(x[0], a3[0]), (a3[0], a2[0]), (a2[0], a4[0]), (a4[0], x[0]), ])

    V = d[1].dense_matrix().right_kernel()
    W = d[2].dense_matrix().transpose().image()
    print(f'{d[1].dimensions() = }, {d[2].dimensions() = }')
    Q = V.quotient(W)
    print(Q)

    v = ket0 - 2 * ket1

    print(f'{v in V = }, {v in W = }')


def join_keep_names(G1, G2):
    G = G1.union(G2)
    G.add_edges([(u, v) for u in G1 for v in G2])
    return G


def demo_fill_cycle_ket0_ketminus():
    G, names = make_graph(2)
    G.plot().save("G.png")
    x, a2, a3, a4, b2, b3, b4 = names
    x1 = '1.x1'

    G_cycle_0 = Graph()
    G_cycle_0.add_edges([(x[0], a3[0]), (a3[0], a2[0]), (a2[0], a4[0]), (a4[0], x[0]), ])
    G_cycle_minus = Graph()
    G_cycle_minus.add_edges(
        [(a2[1], a3[1]), (a3[1], x[1]), (x[1], b3[1]), (b3[1], b2[1]), (b2[1], b4[1]), (b4[1], x1),
         (x1, a4[1]), (a4[1], a2[1]), ])
    G_cycle = join_keep_names(G_cycle_0, G_cycle_minus)
    G_cycle.plot().save("G_cycle.png")
    G_fill = fill_cycle(G, G_cycle, {x1: x[1]})

    for e in G_fill.edges(sort=True):
        print(e[0], e[1])

    S = G_fill.clique_complex()
    print(S.homology())
    C = S.chain_complex(base_ring=QQ)
    d = C.differential()
    print(d)
    L = d[3].transpose() * d[3] + d[4] * d[4].transpose()
    L_ker = L.right_kernel()
    print(L_ker)

    cell2dim = {}
    for i, cell in enumerate(S._n_cells_sorted(3)):
        c, s = sort_and_sign(cell)
        cell2dim[c] = (i, s)

    def bit_cells(i, b):
        return [(x[i], b3[i]), (b3[i], b2[i]), (b2[i], b4[i]), (b4[i], x[i])] if b == 1 \
            else [(x[i], a3[i]), (a3[i], a2[i]), (a2[i], a4[i]), (a4[i], x[i]), ]

    def basis_vector(x, y):
        C0 = bit_cells(0, x)
        C1 = bit_cells(1, y)
        print('LEN', len(C0) * len(C1))
        return cell2vec(cell2dim, [c0 + c1 for c0 in C0 for c1 in C1])

    psi = basis_vector(0, 0) - basis_vector(0, 1)
    phi = basis_vector(0, 0) + basis_vector(0, 1)
    print(f'{psi.inner_product(phi) = }')
    print(f'{basis_vector(0, 0).inner_product(basis_vector(1, 1)) = }')
    print(f'{basis_vector(0, 0).inner_product(basis_vector(0, 1)) = }')

    assert d[3] * phi == 0

    V = d[3].dense_matrix().right_kernel()
    W = d[4].dense_matrix().transpose().image()
    print(f'{d[3].dimensions() = }, {d[4].dimensions() = }')
    Q = V.quotient(W)
    print(Q)

    print(f'{psi in V = }, {psi in W = }')
    print(f'{phi in V = }, {phi in W = }')

    q = Q.quotient_map()
    V2 = span([
        q(basis_vector(0, 0) + basis_vector(0, 1)),
        q(basis_vector(1, 0)),
        q(basis_vector(1, 1)),
    ])
    print(V2.dimension())


def state_00m11():
    G, names = make_graph(2)
    x, a2, a3, a4, b2, b3, b4 = names

    dx = [f'x{i}' for i in range(1, 5)]  # dummy x
    G_cycle = Graph()
    G_cycle.add_edges([(a, b) for a in x for b in dx] + [(a, b) for a, b in zip(dx, dx[1:] + dx[:1])])  # Figure 12

    def bit_cells(i, b):
        return [(x[i], b3[i]), (b3[i], b2[i]), (b2[i], b4[i]), (b4[i], x[i])] if b == 1 \
            else [(x[i], a3[i]), (a3[i], a2[i]), (a2[i], a4[i]), (a4[i], x[i]), ]

    G_cycle_00 = join_keep_names(Graph(bit_cells(0, 0)), Graph(bit_cells(1, 0)))
    loop = [a3[0], a4[1], a4[0], a3[1]]
    G_cycle_00.add_edges([(loop[(i - j) % 4], dx[i]) for i in range(4) for j in range(2)])
    G_cycle_11 = join_keep_names(Graph(bit_cells(0, 1)), Graph(bit_cells(1, 1)))
    loop = [b3[0], b4[1], b4[0], b3[1]]
    G_cycle_11.add_edges([(loop[(i - j) % 4], dx[i]) for i in range(4) for j in range(2)])

    G_cycle = G_cycle.union(G_cycle_00).union(G_cycle_11)
    G_cycle.delete_edge(x[0], x[1])
    G_fill = fill_cycle(G, G_cycle, {xi: x[0] for xi in dx})

    S = G_fill.clique_complex()
    assert str(S.homology()) == '{0: 0, 1: 0, 2: 0, 3: Z x Z x Z, 4: 0}'

    C = S.chain_complex(base_ring=QQ)
    d = C.differential()

    cell2dim = {}
    for i, cell in enumerate(S._n_cells_sorted(3)):
        c, s = sort_and_sign(cell)
        cell2dim[c] = (i, s)

    def basis_vector(x, y):
        C0 = bit_cells(0, x)
        C1 = bit_cells(1, y)
        return cell2vec(cell2dim, [c0 + c1 for c0 in C0 for c1 in C1])

    V = d[3].dense_matrix().right_kernel()
    W = d[4].dense_matrix().transpose().image()
    Q = V.quotient(W)
    assert Q.dimension() == 3

    v = basis_vector(0, 0) - basis_vector(1, 1)
    assert v in V and v in W
    w = basis_vector(0, 0) + basis_vector(1, 1)
    assert w in V and w not in W

    q = Q.quotient_map()
    V2 = span([
        q(basis_vector(0, 0) + basis_vector(1, 1)),
        q(basis_vector(1, 0)),
        q(basis_vector(0, 1)),
    ])
    assert V2.dimension() == 3


def state_01m10():
    G, names = make_graph(2)
    x, a2, a3, a4, b2, b3, b4 = names

    dx = [f'x{i}' for i in range(1, 5)]  # dummy x
    G_cycle = Graph()
    G_cycle.add_edges([(a, b) for a in x for b in dx] + [(a, b) for a, b in zip(dx, dx[1:] + dx[:1])])  # Figure 12

    def bit_cells(i, b):
        return [(x[i], b3[i]), (b3[i], b2[i]), (b2[i], b4[i]), (b4[i], x[i])] if b == 1 \
            else [(x[i], a3[i]), (a3[i], a2[i]), (a2[i], a4[i]), (a4[i], x[i]), ]

    G_cycle_01 = join_keep_names(Graph(bit_cells(0, 0)), Graph(bit_cells(1, 1)))
    loop = [a3[0], a4[1], a4[0], a3[1]]
    G_cycle_01.add_edges([(loop[(i - j) % 4], dx[i]) for i in range(4) for j in range(2)])
    G_cycle_10 = join_keep_names(Graph(bit_cells(0, 1)), Graph(bit_cells(1, 0)))
    loop = [b3[0], b4[1], b4[0], b3[1]]
    G_cycle_10.add_edges([(loop[(i - j) % 4], dx[i]) for i in range(4) for j in range(2)])

    G_cycle = G_cycle.union(G_cycle_01).union(G_cycle_10)
    G_cycle.delete_edge(x[0], x[1])
    G_fill = fill_cycle(G, G_cycle, {xi: x[0] for xi in dx})

    S = G_fill.clique_complex()
    assert str(S.homology()) == '{0: 0, 1: 0, 2: 0, 3: Z x Z x Z, 4: 0}'

    C = S.chain_complex(base_ring=QQ)
    d = C.differential()

    cell2dim = {}
    for i, cell in enumerate(S._n_cells_sorted(3)):
        c, s = sort_and_sign(cell)
        cell2dim[c] = (i, s)

    def basis_vector(x, y):
        C0 = bit_cells(0, x)
        C1 = bit_cells(1, y)
        return cell2vec(cell2dim, [c0 + c1 for c0 in C0 for c1 in C1])

    V = d[3].dense_matrix().right_kernel()
    W = d[4].dense_matrix().transpose().image()
    Q = V.quotient(W)
    assert Q.dimension() == 3

    v = basis_vector(0, 1) - basis_vector(1, 0)
    assert v in V and v in W
    w = basis_vector(0, 1) + basis_vector(1, 0)
    assert w in V and w not in W

    q = Q.quotient_map()
    V2 = span([
        q(basis_vector(0, 1) + basis_vector(1, 0)),
        q(basis_vector(1, 1)),
        q(basis_vector(0, 0)),
    ])
    assert V2.dimension() == 3


def state_00m10m11():
    G, names = make_graph(2)
    x, a2, a3, a4, b2, b3, b4 = names
    x_, a2_, a3_, a4_, b2_, b3_, b4_ = [[f'{v}_' for v in vv] for vv in names]

    dx = [f'x{i}' for i in range(1, 7)]  # dummy x
    G_cycle = Graph()

    def bit_cells(i, b):
        return [(x[i], b3[i]), (b3[i], b2[i]), (b2[i], b4[i]), (b4[i], x[i])] if b == 1 \
            else [(x[i], a3[i]), (a3[i], a2[i]), (a2[i], a4[i]), (a4[i], x[i]), ]

    def bit_cells_(i, b):  # use
        return [(x[i], b3_[i]), (b3_[i], b2_[i]), (b2_[i], b4_[i]), (b4_[i], x[i])] if b == 1 \
            else [(x[i], a3_[i]), (a3_[i], a2_[i]), (a2_[i], a4_[i]), (a4_[i], x[i]), ]

    def make_part(cells0, cells1, dx, loop):
        Gc = join_keep_names(Graph(cells0), Graph(cells1))
        Gc.add_edges([(a, b) for a in x for b in dx] + [(a, b) for a, b in zip(dx, dx[1:] + dx[:1])])  # Figure 12
        Gc.add_edges([(loop[(i - j) % 4], dx[i]) for i in range(4) for j in range(2)])
        return Gc

    G_cycle = G_cycle.union(make_part(bit_cells(0, 0), bit_cells(1, 0), dx[:4], [a3[0], a4[1], a4[0], a3[1]]))
    G_cycle = G_cycle.union(
        make_part(bit_cells(0, 1), bit_cells(1, 1), [dx[1], dx[2], dx[5], dx[4]], [b3[0], b4[1], b4[0], b3[1]]))
    G_cycle = G_cycle.union(
        make_part(bit_cells_(0, 1), bit_cells_(1, 0), [dx[4], dx[5], dx[3], dx[0]], [b3_[0], a4_[1], b4_[0], a3_[1]]))

    G_cycle.delete_edge(x[0], x[1])

    f = {xi: x[0] for xi in dx}
    for vv, vv_ in [(a2, a2_), (a3, a3_), (a4, a4_), (b2, b2_), (b3, b3_), (b4, b4_)]:
        for v, v_ in zip(vv, vv_):
            f[v_] = v
    G_fill = fill_cycle(G, G_cycle, f)

    S = G_fill.clique_complex()
    assert str(S.homology()) == '{0: 0, 1: 0, 2: 0, 3: Z x Z x Z, 4: 0, 5: 0}'

    C = S.chain_complex(base_ring=QQ)
    d = C.differential()

    cell2dim = {}
    for i, cell in enumerate(S._n_cells_sorted(3)):
        c, s = sort_and_sign(cell)
        cell2dim[c] = (i, s)

    def basis_vector(x, y):
        C0 = bit_cells(0, x)
        C1 = bit_cells(1, y)
        return cell2vec(cell2dim, [c0 + c1 for c0 in C0 for c1 in C1])

    V = d[3].dense_matrix().right_kernel()
    W = d[4].dense_matrix().transpose().image()
    Q = V.quotient(W)
    assert Q.dimension() == 3

    v = basis_vector(0, 0) - basis_vector(1, 0) - basis_vector(1, 1)
    assert v in V and v in W

    q = Q.quotient_map()
    V2 = span([
        q(basis_vector(0, 0) + basis_vector(1, 1)),
        q(basis_vector(1, 0) - basis_vector(1, 1)),
        q(basis_vector(0, 1)),
    ])
    assert V2.dimension() == 3


def state_01m10p11():
    G, names = make_graph(2)
    x, a2, a3, a4, b2, b3, b4 = names
    x_, a2_, a3_, a4_, b2_, b3_, b4_ = [[f'{v}_' for v in vv] for vv in names]

    dx = [f'x{i}' for i in range(1, 7)]  # dummy x
    G_cycle = Graph()

    def bit_cells(i, b):
        return [(x[i], b3[i]), (b3[i], b2[i]), (b2[i], b4[i]), (b4[i], x[i])] if b == 1 \
            else [(x[i], a3[i]), (a3[i], a2[i]), (a2[i], a4[i]), (a4[i], x[i]), ]

    def bit_cells_(i, b):  # use
        return [(x[i], b3_[i]), (b3_[i], b2_[i]), (b2_[i], b4_[i]), (b4_[i], x[i])] if b == 1 \
            else [(x[i], a3_[i]), (a3_[i], a2_[i]), (a2_[i], a4_[i]), (a4_[i], x[i]), ]

    def make_part(cells0, cells1, dx, loop):
        Gc = join_keep_names(Graph(cells0), Graph(cells1))
        Gc.add_edges([(a, b) for a in x for b in dx] + [(a, b) for a, b in zip(dx, dx[1:] + dx[:1])])  # Figure 12
        Gc.add_edges([(loop[(i - j) % 4], dx[i]) for i in range(4) for j in range(2)])
        return Gc

    G_cycle = G_cycle.union(make_part(bit_cells(0, 0), bit_cells(1, 1), dx[:4], [a3[0], b4[1], a4[0], b3[1]]))
    G_cycle = G_cycle.union(
        make_part(bit_cells(0, 1), bit_cells(1, 0), [dx[1], dx[2], dx[5], dx[4]], [b3[0], a4[1], b4[0], a3[1]]))
    G_cycle = G_cycle.union(
        make_part(bit_cells_(0, 1), bit_cells_(1, 1), [dx[4], dx[0], dx[3], dx[5]], [b3_[0], b4_[1], b4_[0], b3_[1]]))

    G_cycle.delete_edge(x[0], x[1])

    f = {xi: x[0] for xi in dx}
    for vv, vv_ in [(a2, a2_), (a3, a3_), (a4, a4_), (b2, b2_), (b3, b3_), (b4, b4_)]:
        for v, v_ in zip(vv, vv_):
            f[v_] = v
    G_fill = fill_cycle(G, G_cycle, f)

    S = G_fill.clique_complex()
    assert str(S.homology()) == '{0: 0, 1: 0, 2: 0, 3: Z x Z x Z, 4: 0, 5: 0}'

    C = S.chain_complex(base_ring=QQ)
    d = C.differential()

    cell2dim = {}
    for i, cell in enumerate(S._n_cells_sorted(3)):
        c, s = sort_and_sign(cell)
        cell2dim[c] = (i, s)

    def basis_vector(x, y):
        C0 = bit_cells(0, x)
        C1 = bit_cells(1, y)
        return cell2vec(cell2dim, [c0 + c1 for c0 in C0 for c1 in C1])

    V = d[3].dense_matrix().right_kernel()
    W = d[4].dense_matrix().transpose().image()
    Q = V.quotient(W)
    assert Q.dimension() == 3

    v = basis_vector(0, 1) - basis_vector(1, 0) + basis_vector(1, 1)
    assert v in V and v in W

    q = Q.quotient_map()
    V2 = span([
        q(basis_vector(0, 1) + basis_vector(1, 0)),
        q(basis_vector(1, 0) + basis_vector(1, 1)),
        q(basis_vector(0, 0)),
    ])
    assert V2.dimension() == 3


def state_000():
    G, names = make_graph(3)
    x, a2, a3, a4, b2, b3, b4 = names

    def bit_cells(i, b):
        return [(x[i], b3[i]), (b3[i], b2[i]), (b2[i], b4[i]), (b4[i], x[i])] if b == 1 \
            else [(x[i], a3[i]), (a3[i], a2[i]), (a2[i], a4[i]), (a4[i], x[i]), ]

    G_cycle = join_keep_names(Graph(bit_cells(0, 0)), Graph(bit_cells(1, 0)))
    G_cycle = join_keep_names(G_cycle, Graph(bit_cells(2, 0)))
    G_fill = fill_cycle(G, G_cycle)

    S = G_fill.clique_complex()
    assert str(S.homology()) == '{0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: Z^7, 6: 0}'
    C = S.chain_complex(base_ring=QQ)
    d = C.differential()

    cell2dim = {}
    for i, cell in enumerate(S._n_cells_sorted(5)):
        c, s = sort_and_sign(cell)
        cell2dim[c] = (i, s)

    def basis_vector(x, y, z):
        C0 = bit_cells(0, x)
        C1 = bit_cells(1, y)
        C2 = bit_cells(2, z)
        return cell2vec(cell2dim, [c0 + c1 + c2 for c0 in C0 for c1 in C1 for c2 in C2])

    V = d[5].dense_matrix().right_kernel()
    W = d[6].dense_matrix().transpose().image()
    Q = V.quotient(W)
    assert Q.dimension() == 7

    v = basis_vector(0, 0, 0)
    assert v in V and v in W

    q = Q.quotient_map()
    V2 = span([
        q(basis_vector(0, 0, 1)),
        q(basis_vector(0, 1, 0)),
        q(basis_vector(0, 1, 1)),
        q(basis_vector(1, 0, 0)),
        q(basis_vector(1, 0, 1)),
        q(basis_vector(1, 1, 0)),
        q(basis_vector(1, 1, 1)),
    ])
    assert V2.dimension() == 7


state_0m1()
state_00m11()
state_01m10()
state_00m10m11()
state_01m10p11()
state_000()
