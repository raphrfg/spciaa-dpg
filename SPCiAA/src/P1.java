
/**
 * @(#)P1.java
 *
 *
 * @author 
 * @version 1.00 2015/2/23
 */
import java.util.concurrent.CyclicBarrier;

class P1 extends Production {
    P1(Vertex Vert, CyclicBarrier Barrier) {
        super(Vert, Barrier);
    }

    Vertex apply(Vertex S) {
        System.out.println("p1");
        Vertex T1 = new Vertex(null, null, S, "T");
        Vertex T2 = new Vertex(null, null, S, "T");
        S.set_left(T1);
        S.set_right(T2);
        S.set_label("root");
        return S;
    }
}
