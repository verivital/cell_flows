% show a plot of the routing tree of next pointers for each color
for tt = 1 : length(targets)
    figure;
    hold on;
    treeplot(tree(:,tt)','', colors(mod(tt, length(colors)) + 1));
    count = size(tree(:,tt)',2);
    [x,y] = treelayout(tree(:,tt)');
    x = x';
    y = y';
    name1 = cellstr(num2str((1:count)'));
    text(x(:,1), y(:,1), name1, 'VerticalAlignment','bottom','HorizontalAlignment','right')
end