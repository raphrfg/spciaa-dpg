
public class Main {

    public static void main(String[] args) throws Exception {
        Executor s = new Executor();
        
        int k = 7;
        double dt = 0.01;
        int steps = 100;
        s.run(k, dt, steps);
    }

}