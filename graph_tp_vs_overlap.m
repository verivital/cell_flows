tp = [0.1456	0.0744	0.0296	0.0184	0.0136	0.0108	0.0088	0.0076	0.0076]';

path_length = 8;
pc = [0:path_length]/path_length;

plot(pc, tp, 'Color', [205 92 92]/255, 'LineWidth', 2);
xlabel('fraction of overlapping paths');
ylabel('throughput');
ylim([0 0.15]);


figure;
hold on;
path_length = [1	2	3	4	5	6	7	8	9	10	11	12	13	14	15];
tp_square = [0.2492	0.1412	0.1052	0.0812	0.076	0.0756	0.0748	0.074	0.0732	0.0724	0.0716	0.0708	0.07	0.0692	0.0684];
tp_eqtri = [0.2496	0.1108	0.0828	0.066	0.0656	0.0656	0.0652	0.0648	0.0648	0.0644	0.064	0.064	0.0636	0.0632	0.0632];

plot(path_length, tp_square, '-s', 'Color', [70 130 180]/255, 'LineWidth', 2);
plot(path_length, tp_eqtri, '-v', 'Color', [205 92 92]/255, 'LineWidth', 2);
xlabel('path length');
ylabel('throughput');
legend('square','triangular');
%ylim([0 0.15]);


figure;
hold on;
color_overlap = [0:8];

tp_color = [0.1412	0.0316	0.0316	0.0316	0.0316	0.0316	0.0316	0.0316	0.032];
plot(color_overlap, tp_color, 'Color', [205 92 92]/255, 'LineWidth', 2);
xlabel('number of overlapping colors');
ylabel('throughput');
%ylim([0 0.15]);