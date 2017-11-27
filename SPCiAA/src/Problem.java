import java.util.concurrent.CyclicBarrier;

public interface Problem {

    Production makeA1(Vertex Vert, double h, double dt, double t, CyclicBarrier Barrier);

    Production makeA(Vertex Vert, double h, double dt, double t, CyclicBarrier Barrier);

    Production makeAN(Vertex Vert, double h, double dt, double t, CyclicBarrier Barrier);
}
