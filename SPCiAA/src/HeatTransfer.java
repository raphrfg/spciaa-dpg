import java.util.concurrent.CyclicBarrier;

public class HeatTransfer implements Problem {

    @Override
    public Production makeA1(Vertex v, double h, double dt, double t, CyclicBarrier barrier) {
        return new ABase(v, h, dt, t, barrier) {
            @Override
            Vertex apply(Vertex v) {
                double a = dt / (h * h);
                v.m_a[1][1] = 1.0;
                v.m_a[1][2] = 0.0;
                v.m_b[1] = 0;

                v.m_a[2][1] = -a;
                v.m_a[2][2] = 0.5 + a;
                v.m_b[2] = 0.5 * v.m_x[2];
                return v;
            }
        };
    }

    @Override
    public Production makeA(Vertex v, double h, double dt, double t, CyclicBarrier barrier) {
        return new ABase(v, h, dt, t, barrier) {
            @Override
            Vertex apply(Vertex v) {
                double a = dt / (h * h);
                v.m_a[1][1] = 0.5 + a;
                v.m_a[1][2] = -a;
                v.m_b[1] = 0.5 * v.m_x[1];

                v.m_a[2][1] = -a;
                v.m_a[2][2] = 0.5 + a;
                v.m_b[2] = 0.5 * v.m_x[2];
                return v;
            }
        };
    }

    @Override
    public Production makeAN(Vertex v, double h, double dt, double t, CyclicBarrier barrier) {
        return new ABase(v, h, dt, t, barrier) {
            @Override
            Vertex apply(Vertex v) {
                double a = dt / (h * h);
                v.m_a[1][1] = 0.5 + a;
                v.m_a[1][2] = -a;
                v.m_b[1] = 0.5 * v.m_x[1];

                v.m_a[2][1] = 0.0;
                v.m_a[2][2] = 1.0;
                v.m_b[2] = 0;
                return v;
            }
        };
    }

}
