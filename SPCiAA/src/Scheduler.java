import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;

public class Scheduler {

    private CyclicBarrier barrier;
    private List<Production> queue = new ArrayList<>();

    public void beginStage(int productions) {
        barrier = new CyclicBarrier(productions + 1);
        queue.clear();
    }

    public void add(Production prod) {
        queue.add(prod);
    }

    public void executeStage() throws InterruptedException, BrokenBarrierException {
        for (final Production prod : queue) {
            prod.start();
        }
        barrier.await();
    }

    public void executeOne(Function<CyclicBarrier, Production> makeProd)
            throws InterruptedException, BrokenBarrierException {
        beginStage(1);
        add(makeProd.apply(barrier));
        executeStage();
    }

    public void forEachVertex(Collection<Vertex> vs, BiFunction<Vertex, CyclicBarrier, Production> makeProd)
            throws InterruptedException, BrokenBarrierException {
        beginStage(vs.size());
        for (Vertex v : vs) {
            add(makeProd.apply(v, barrier));
        }
        executeStage();
    }

    public CyclicBarrier barrier() {
        return barrier;
    }

    public List<Production> productions() {
        return Collections.unmodifiableList(queue);
    }

}
