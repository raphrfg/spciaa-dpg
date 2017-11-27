import java.util.concurrent.CyclicBarrier;

class A2 extends Production {
    A2(Vertex Vert, CyclicBarrier Barrier) {
        super(Vert, Barrier);
    }

    Vertex apply(Vertex T) {
        System.out.println("A2");
        T.m_a[0][0] = T.m_left.m_a[2][2] + T.m_right.m_a[1][1];
        T.m_a[1][0] = T.m_left.m_a[1][2];
        T.m_a[2][0] = T.m_right.m_a[2][1];
        T.m_a[0][1] = T.m_left.m_a[2][1];
        T.m_a[1][1] = T.m_left.m_a[1][1];
        T.m_a[2][1] = 0.0;
        T.m_a[0][2] = T.m_right.m_a[1][2];
        T.m_a[1][2] = 0.0;
        T.m_a[2][2] = T.m_right.m_a[2][2];
        T.m_b[0] = T.m_left.m_b[2] + T.m_right.m_b[1];
        T.m_b[1] = T.m_left.m_b[1];
        T.m_b[2] = T.m_right.m_b[2];
        return T;
    }
}
