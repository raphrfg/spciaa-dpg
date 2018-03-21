import java.util.concurrent.CyclicBarrier;

public class AdvectionExplicit implements Problem {
    
    private static final double beta = 1;
    private double F(double t) {
        return 2*Math.PI*Math.cos(2 * Math.PI * t);
    }

    @Override
    public Production makeA1(Vertex v, double h, double dt, double t, CyclicBarrier barrier) {
        return new ABase(v, h, dt, t, barrier) {
            @Override
            Vertex apply(Vertex v) {

                double a = beta* dt / (2 * h);
                v.m_a[1][1] = 1.0;
                v.m_a[1][2] = 0.0;
                v.m_b[1] = 0 ;

                v.m_a[2][1] = 0;
                v.m_a[2][2] = 0.5;
                v.m_b[2] = v.m_x[2]+a*v.m_x[1]+dt*F(t);

                return v;
            }
        };
    }

    @Override
    public Production makeA(Vertex v, double h, double dt, double t, CyclicBarrier barrier) {
        return new ABase(v, h, dt, t, barrier) {
            @Override
            Vertex apply(Vertex v) {
                double a = beta* dt / (2 * h);
                v.m_a[1][1] = 0.5;
                v.m_a[1][2] = 0.0;
                v.m_b[1] = -a*v.m_x[2];

                v.m_a[2][1] = 0;
                v.m_a[2][2] = 0.5;
                v.m_b[2] = v.m_x[2] + a * v.m_x[1] + dt * F(t);
                return v;
            }
        };
    }

    @Override
    public Production makeAN(Vertex v, double h, double dt, double t, CyclicBarrier barrier) {
        return new ABase(v, h, dt, t, barrier) {
            @Override
            Vertex apply(Vertex v) {
                double a = beta* dt / (2 * h);
                v.m_a[1][1] = 0.5;
                v.m_a[1][2] = 0.0;
                v.m_b[1] = -a*v.m_x[2];
                v.m_a[2][1] = 0;
                v.m_a[2][2] = 1;
                v.m_b[2] = (1-2*a)*v.m_x[2] + 2*a * v.m_x[1] + dt * F(t);
                return v;
            }
        };
    }

}
