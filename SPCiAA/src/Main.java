
public class Main {

    public static void main(String[] args) throws Exception {
        Executor s = new Executor();
        
        int k = 8;
        double dt = 0.01;
        int steps = 1000;
        s.run(k, dt, steps);
    }

}