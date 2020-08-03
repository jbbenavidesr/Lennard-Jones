/**
 * Defines useful functions for animationg with gnuplot.
 * @param animation: If 2, plot without saving. If 3, saves in a gif
 */

void start_animation(int animation)
{
    if (animation > 1)
    {
        if (animation == 3)
        {
            std::cout << "set terminal gif animate\n";
            std::cout << "set output 'movie01.gif'\n";
        }
        std::cout << "unset key\n";
        std::cout << "set xrange[0:" << Lx << "]\n";
        std::cout << "set yrange[0:" << Lz << "]\n";
        std::cout << "set size ratio -1\n";
        std::cout << "set parametric\n";
        std::cout << "set trange [0:7]\n";
        std::cout << "set isosamples 12" << std::endl;
    }
}
void begin_frame(int animation)
{
    if (animation > 1)
    {
        std::cout << "plot 0,0 ";
        //std::cout << " , " << Lx/7 << "*t,0";             // down wall
        //std::cout << " , " << Lx/7 << "*t," << Ly;        // up wall
        //std::cout << " , 0," << Ly/7 << "*t";             // left wall
        //std::cout << " , " << Lx << "," << Ly/7 << "*t";  // right wall
    }
}
void end_frame(int animation)
{
    if (animation > 1)
        std::cout << std::endl;
}

void plot_single_point(int animation, double x, double y)
{
    if (animation > 1)
        std::cout << " , '< echo " << x << ' ' << y << "' w points";
}
