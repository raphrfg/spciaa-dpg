import java.util.concurrent.CyclicBarrier;

class BS extends Production {
    BS(Vertex Vert, CyclicBarrier Barrier) {
        super(Vert, Barrier);
    }

    Vertex apply(Vertex T) {
        System.out.println("BS");
        T.m_left.m_x[0] = T.m_x[0];
        T.m_left.m_x[2] = T.m_x[1];
        T.m_left.m_x[1] = (T.m_left.m_b[0] - T.m_left.m_a[0][1] * T.m_left.m_x[0]
                - T.m_left.m_a[0][2] * T.m_left.m_x[2]) / T.m_left.m_a[0][0];
        T.m_right.m_x[0] = T.m_x[1];
        T.m_right.m_x[2] = T.m_x[2];
        T.m_right.m_x[1] = (T.m_right.m_b[0] - T.m_right.m_a[0][1] * T.m_right.m_x[0]
                - T.m_right.m_a[0][2] * T.m_right.m_x[2]) / T.m_right.m_a[0][0];
        return T;
    }
}