import java.util.concurrent.CyclicBarrier;

class BSA extends Production {
    BSA(Vertex Vert, CyclicBarrier Barrier) {
        super(Vert, Barrier);
    }

    Vertex apply(Vertex T) {
        System.out.println("BSA");
        T.m_left.m_x[1] = T.m_x[0];
        T.m_left.m_x[2] = T.m_x[1];

        T.m_right.m_x[1] = T.m_x[1];
        T.m_right.m_x[2] = T.m_x[2];
        return T;
    }
}