import java.util.concurrent.CyclicBarrier;

public abstract class ABase extends Production {

    protected double h;
    protected double dt;
    protected double t;

    ABase(Vertex Vert, double h, double dt, double t, CyclicBarrier barrier) {
        super(Vert, barrier);
        this.h = h;
        this.dt = dt;
        this.t = t;
    }

}
