import java.util.concurrent.CyclicBarrier;

class E2 extends Production {
    E2(Vertex Vert, CyclicBarrier Barrier) {
        super(Vert, Barrier);
    }

    Vertex apply(Vertex T) {
        System.out.println("E2");
        T.m_b[0] /= T.m_a[0][0];
        T.m_a[0][2] /= T.m_a[0][0];
        T.m_a[0][1] /= T.m_a[0][0];
        T.m_a[0][0] /= T.m_a[0][0];
        T.m_b[1] -= T.m_b[0] * T.m_a[1][0];
        T.m_a[1][2] -= T.m_a[0][2] * T.m_a[1][0];
        T.m_a[1][1] -= T.m_a[0][1] * T.m_a[1][0];
        T.m_a[1][0] -= T.m_a[0][0] * T.m_a[1][0];
        T.m_b[2] -= T.m_b[0] * T.m_a[2][0];
        T.m_a[2][2] -= T.m_a[0][2] * T.m_a[2][0];
        T.m_a[2][1] -= T.m_a[0][1] * T.m_a[2][0];
        T.m_a[2][0] -= T.m_a[0][0] * T.m_a[2][0];
        return T;
    }
}