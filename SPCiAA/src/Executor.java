import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toSet;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

class Executor {

    final Scheduler sched = new Scheduler();
    final Vertex root = new Vertex(null, null, null, "S");
    
    final Problem problem = new AdvectionImplicit();

    private Set<Vertex> parentsOf(Collection<Vertex> vs) {
        return vs.stream().map(v -> v.m_parent).collect(toSet());
    }

    public void run(int k, double dt, int steps) throws Exception {
        final List<Vertex> leafs = buildTree(k);

        Function<Double, Double> init = x -> {
//            double r = 5 * abs(x - 0.25);
//            return r > 1 ? 0 : (r + 1) * (r + 1) * (r - 1) * (r - 1);
            return Math.sin(2 * Math.PI *x);
        };
        setInitState(leafs, init);

        plotSolution(leafs);

        for (int i = 0; i < steps; ++i) {
            solverStep(k, dt, i * dt, leafs);
            plotSolution(leafs);
        }
    }

    private void solverStep(int k, double dt, double t, final List<Vertex> leafs)
            throws InterruptedException, BrokenBarrierException {
        computeLeafMatrices(k, dt, t, leafs);

        // going up
        Set<Vertex> level = parentsOf(leafs);
        while (level.size() > 1) {
            sched.forEachVertex(level, A2::new);
            sched.forEachVertex(level, E2::new);
            level = parentsOf(level);
        }
        sched.executeOne(b -> new Aroot(root, b));
        sched.executeOne(b -> new Eroot(root, b));

        backwardSub(k);
    }

    private void setInitState(List<Vertex> leafs, Function<Double, Double> init) {
        double h = 1.0 / (leafs.size() + 1);
        for (int i = 0; i < leafs.size(); ++i) {
            double x1 = i * h;
            double x2 = (i + 1) * h;
            double x3 = (i + 2) * h;
            leafs.get(i).m_x[1] = init.apply(x1);
            leafs.get(i).m_x[2] = init.apply(x2);
        }
    }

    private List<Vertex> buildTree(int k) throws InterruptedException, BrokenBarrierException {
        sched.executeOne(b -> new P1(root, b));

        // internal nodes
        List<Vertex> previousLevel = Arrays.asList(root);
        for (int i = 1; i < k - 1; ++i) {
            sched.beginStage(2 * previousLevel.size());
            for (Vertex v : previousLevel) {
                sched.add(new P2(v.m_left, sched.barrier()));
                sched.add(new P2(v.m_right, sched.barrier()));
            }
            sched.executeStage();
            previousLevel = sched.productions().stream().map(p -> p.m_vertex).collect(toList());
        }

        // leafs
        sched.beginStage(2 * previousLevel.size());
        for (Vertex v : previousLevel) {
            sched.add(new P3(v.m_left, sched.barrier()));
            sched.add(new P3(v.m_right, sched.barrier()));
        }
        sched.executeStage();
        return sched.productions().stream().map(p -> p.m_vertex).collect(toList());
    }

    private void computeLeafMatrices(int k, double dt, double t, final List<Vertex> leafs)
            throws InterruptedException, BrokenBarrierException {
        double h = Math.pow(2, -(k - 1));

        sched.beginStage(leafs.size());
        
        sched.add(problem.makeA1(leafs.get(0), h, dt, t, sched.barrier()));
        for (int i = 1; i < leafs.size() - 1; ++ i) {
            sched.add(problem.makeA(leafs.get(i), h, dt, t, sched.barrier()));
        }
        sched.add(problem.makeAN(leafs.get(leafs.size() - 1), h, dt, t, sched.barrier()));
        
        sched.executeStage();
    }
    
    private void backwardSub(int k) throws InterruptedException, BrokenBarrierException {
        sched.executeOne(b -> new BS(root, b));

        List<Vertex> previousLevel = Arrays.asList(root);
        for (int i = 1; i < k - 2; ++i) {
            sched.beginStage(2 * previousLevel.size());
            for (Vertex v : previousLevel) {
                sched.add(new BS(v.m_left, sched.barrier()));
                sched.add(new BS(v.m_right, sched.barrier()));
            }
            sched.executeStage();
            previousLevel = sched.productions().stream().map(p -> p.m_vertex).collect(toList());
        }

        sched.beginStage(2 * previousLevel.size());
        for (Vertex v : previousLevel) {
            sched.add(new BSA(v.m_left, sched.barrier()));
            sched.add(new BSA(v.m_right, sched.barrier()));
        }
        sched.executeStage();
    }

    private void plotSolution(final List<Vertex> leafs) throws InterruptedException {
        int delay = 300;
        double[] initState = extractSolution(leafs);
        ResultPrinter.printResult(initState);
//        ResultPrinter.plot.setFixedBounds(1, -0.2, 1.8);
        ResultPrinter.plot.setFixedBounds(1, -1.5, 1.5);
        TimeUnit.MILLISECONDS.sleep(delay);
    }

    private double[] extractSolution(final List<Vertex> leafs) {
        double[] solution = new double[leafs.size() + 1];
        solution[0] = leafs.get(0).m_x[1];
        for (int i = 0; i < leafs.size(); ++i) {
            solution[i + 1] = leafs.get(i).m_x[2];
        }
        return solution;
    }
}
