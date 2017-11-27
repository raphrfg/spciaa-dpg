
import java.awt.Dimension;
import java.lang.reflect.InvocationTargetException;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.math.plot.Plot2DPanel;

public class ResultPrinter {

    public static final Plot2DPanel plot = new Plot2DPanel();
    private static final JFrame frame = new JFrame("Plot panel");

    private static final Dimension SIZE = new Dimension(1000, 650);

    static {
        plot.setPreferredSize(SIZE);
        frame.add(plot);
        frame.pack();
        frame.setVisible(true);
    }

    public static void printResult(double[] result) {
        try {
            SwingUtilities.invokeAndWait(() -> {
                double[] x = new double[result.length];
                for (int i = 0; i < result.length; ++i) {
                    x[i] = i / (double) (result.length - 1);
                }
                plot.removeAllPlots();
                plot.addLinePlot("my plot", x, result);
            });
        } catch (Exception e) {
            System.err.println("Plotting issues");
            e.printStackTrace(System.err);
        }
    }
}
