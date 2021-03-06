import static java.lang.Math.PI;

import java.util.concurrent.CyclicBarrier;

public class AdvectionDPG implements Problem {
    
    private static final double beta = 1;
    
    private double F(double t) {
        return 2 * PI * Math.cos(2 * PI * t);
    }

    @Override
    public Production makeA1(Vertex v, double h, double dt, double t, CyclicBarrier barrier) {
        return new ABase(v, h, dt, t, barrier) {
            @Override
            Vertex apply(Vertex T) {
                return T;
            }
        };
    }

    @Override
    public Production makeA(Vertex v, double h, double dt, double t, CyclicBarrier barrier) {
        return new ABase(v, h, dt, t, barrier) {
            @Override
            Vertex apply(Vertex T) {
                return T;
            }
        };
    }

    @Override
    public Production makeAN(Vertex v, double h, double dt, double t, CyclicBarrier barrier) {
        return new ABase(v, h, dt, t, barrier) {
            @Override
            Vertex apply(Vertex T) {
                return T;
            }
        };
    }

}
