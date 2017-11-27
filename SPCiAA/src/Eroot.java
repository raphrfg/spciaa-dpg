import java.util.concurrent.CyclicBarrier;

class Eroot extends Production {
    Eroot(Vertex Vert, CyclicBarrier Count) {
        super(Vert, Count);
    }

    Vertex apply(Vertex T) {
        System.out.println("Eroot");

        // divide first row by diagonal
        T.m_b[1] /= T.m_a[1][1];
        T.m_a[1][2] /= T.m_a[1][1];
        T.m_a[1][0] /= T.m_a[1][1];
        T.m_a[1][1] /= T.m_a[1][1];
        // 2nd=2nd-1st*diag
        T.m_b[0] -= T.m_b[1] * T.m_a[0][1];
        T.m_a[0][2] -= T.m_a[1][2] * T.m_a[0][1];
        T.m_a[0][0] -= T.m_a[1][0] * T.m_a[0][1];
        T.m_a[0][1] -= T.m_a[1][1] * T.m_a[0][1];
        // 3rd=3rd-1st*diag
        T.m_b[2] -= T.m_b[1] * T.m_a[2][1];
        T.m_a[2][2] -= T.m_a[1][2] * T.m_a[2][1];
        T.m_a[2][0] -= T.m_a[1][0] * T.m_a[2][1];
        T.m_a[2][1] -= T.m_a[1][1] * T.m_a[2][1];

        // divide second row by diagonal
        T.m_b[0] /= T.m_a[0][0];
        T.m_a[0][2] /= T.m_a[0][0];
        T.m_a[0][0] /= T.m_a[0][0];

        // 3rd=3rd-2nd*diag
        T.m_b[2] -= T.m_b[0] * T.m_a[2][0];
        T.m_a[2][2] -= T.m_a[0][2] * T.m_a[2][0];
        T.m_a[2][0] -= T.m_a[0][0] * T.m_a[2][0];

        // divide third row by diagonal
        T.m_b[2] /= T.m_a[2][2];
        T.m_a[2][2] /= T.m_a[2][2];

        // b.s.
        T.m_x[2] = T.m_b[2] / T.m_a[2][2];
        T.m_x[1] = (T.m_b[0] - T.m_a[0][2] * T.m_x[2]) / T.m_a[0][0];
        T.m_x[0] = (T.m_b[1] - T.m_a[1][0] * T.m_x[1] - T.m_a[1][2] * T.m_x[2]) / T.m_a[1][1];

        return T;
    }
}