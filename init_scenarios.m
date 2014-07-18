OPT_TESSELATION_TRI = 0;
OPT_TESSELATION_SQUARE = 1;
OPT_TESSELATION_SQUARE_LINE = 2;
OPT_TESSELATION_RECTANGULAR = 3;
OPT_TESSELATION_SNUB_SQUARE_TILING = 4;
OPT_TESSELATION_ELONGATED_TRIANGULAR_TILING = 5;
OPT_TESSELATION_PARALLELOGRAM = 6; % parallelogram works, but vectors are different: move entities with a vector in the same direction as the opposite parallel side (e.g., /_/ going up would move with the vector /)
OPT_TESSELATION_HEXAGON = 7;
OPT_TESSELATION_EQTRI_LINE = 8; % line of equilateral triangles
OPT_TESSELATION_RTRI_LINE = 9; % line of right triangles
OPT_TESSELATION_HEXAGON_DIAMOND = 10;

opt_single = 21; % 0 = batch mode, run all simulations, else value is the number of a simulation below
opt_run = 1; % 0 = batch mode, run all simulations, else value is the number of a run of a simulation below

%MAX_ROUNDS = N + 1;    % needs to be at least N (for routing to have necessarily stabilized), change as necessary
%MAX_ROUNDS = inf;       % run forever

simulations(1).environment = OPT_TESSELATION_SQUARE_LINE;
simulations(1).maxRounds = 2500;
simulations(1).run(1).sources = [1]; % 0 overlap
simulations(1).run(1).targets = [3];

simulations(1).run(2).sources = [1 1]; % 1 overlap
simulations(1).run(2).targets = [3 3];

simulations(1).run(3).sources = [1 1 1]; % 2 overlap
simulations(1).run(3).targets = [3 3 3];

simulations(1).run(4).sources = [1 1 1 1]; % 3 overlap
simulations(1).run(4).targets = [3 3 3 3];

simulations(1).run(5).sources = [1 1 1 1 1]; % 4 overlap
simulations(1).run(5).targets = [3 3 3 3 3];

simulations(1).run(6).sources = [1 1 1 1 1 1]; % 5 overlap
simulations(1).run(6).targets = [3 3 3 3 3 3];

simulations(1).run(7).sources = [1 1 1 1 1 1 1]; % 6 overlap
simulations(1).run(7).targets = [3 3 3 3 3 3 3];

simulations(1).run(8).sources = [1 1 1 1 1 1 1 1]; % 7 overlap
simulations(1).run(8).targets = [3 3 3 3 3 3 3 3];

simulations(1).run(9).sources = [1 1 1 1 1 1 1 1 1]; % 8 overlap
simulations(1).run(9).targets = [3 3 3 3 3 3 3 3 3];


simulations(2).environment = OPT_TESSELATION_SQUARE_LINE;
simulations(2).maxRounds = 2500;
simulations(2).run(1).sources = [8 9]; % square left-right throughput, 0 overlap  (no turns)
simulations(2).run(1).targets = [1 16];

simulations(2).run(2).sources = [8 8]; % square left-right throughput, 1 overlap  (no turns)
simulations(2).run(2).targets = [1 15]; 

simulations(2).run(3).sources = [8 7]; % square left-right throughput, 2 overlap  (no turns)
simulations(2).run(3).targets = [1 14]; 

simulations(2).run(4).sources = [8 6]; % square left-right throughput, 3 overlap  (no turns)
simulations(2).run(4).targets = [1 13];

simulations(2).run(5).sources = [8 5]; % square left-right throughput, 4 overlap  (no turns)
simulations(2).run(5).targets = [1 12];

simulations(2).run(6).sources = [8 4]; % square left-right throughput, 5 overlap  (no turns)
simulations(2).run(6).targets = [1 11];

simulations(2).run(7).sources = [8 3]; % square left-right throughput, 6 overlap  (no turns)
simulations(2).run(7).targets = [1 10];

simulations(2).run(8).sources = [8 2]; % square left-right throughput, 7 overlap  (no turns)
simulations(2).run(8).targets = [1 9];

simulations(2).run(9).sources = [8 1]; % square left-right throughput, 8 overlap  (no turns)
simulations(2).run(9).targets = [1 8];

simulations(3).environment = OPT_TESSELATION_SQUARE_LINE;
simulations(3).maxRounds = 2500;
start = 0;
for tmp = start + 1 : 16
    simulations(3).run(tmp - start).sources = [1];
    simulations(3).run(tmp - start).targets = [tmp];
end

simulations(4).environment = OPT_TESSELATION_EQTRI_LINE;
simulations(4).maxRounds = 500;
start = 0;
for tmp = start + 1 : 16
    simulations(4).run(tmp - start).sources = [1];
    simulations(4).run(tmp - start).targets = [tmp];
end

simulations(5).environment = OPT_TESSELATION_RTRI_LINE;
simulations(5).maxRounds = 500;
start = 0;
for tmp = start + 1 : 16
    simulations(5).run(tmp - start).sources = [1];
    simulations(5).run(tmp - start).targets = [tmp];
end



simulations(7).environment = OPT_TESSELATION_SQUARE;
simulations(7).maxRounds = 2500;
simulations(7).run(1).sources = [3 18]; % square: overlapping
simulations(7).run(1).targets = [23 8];

simulations(8).environment = OPT_TESSELATION_SQUARE;
simulations(8).maxRounds = 2500;
simulations(8).run(1).sources = [3 11]; % square left-right throughput, 2 overlap (has turns)
simulations(8).run(1).targets = [23 19]; 
simulations(8).run(1).failed = [14 9 7];

simulations(9).environment = OPT_TESSELATION_SQUARE;
simulations(9).maxRounds = 2500;
simulations(9).run(1).sources = [3 7]; % square left-right throughput, 3 overlap  (has turns)
simulations(9).run(1).targets = [23 19]; 
simulations(9).run(1).failed = [12 17 9 14];


simulations(10).environment = OPT_TESSELATION_SQUARE;
simulations(10).maxRounds = 2500;
simulations(10).run(1).sources = [3 18 25 21 6]; % square: overlapping, multi-lock disjoint case
simulations(10).run(1).targets = [23 8 5 1 16];

simulations(11).environment = OPT_TESSELATION_SQUARE;
simulations(11).maxRounds = 2500;
simulations(11).run(1).sources = [4 20 22 6]; % square corners
simulations(11).run(1).targets = [25 21 1 5];
simulations(11).run(1).failed = [7 8 9 12 13 14 17 18 19];

simulations(11).run(2).sources = [4 20 22 6 25 21 1 5]; % square corners, 2 copies
simulations(11).run(2).targets = [25 21 1 5 4 20 22 6];
simulations(11).run(3).failed = [7 8 9 12 13 14 17 18 19];

simulations(11).run(3).sources = [4 20 22 6 4 20 22 6 4 20 22 6]; % square corners, 3 copies
simulations(11).run(3).targets = [25 21 1 5 25 21 1 5 25 21 1 5];
simulations(11).run(3).failed = [7 8 9 12 13 14 17 18 19];


simulations(12).environment = OPT_TESSELATION_SQUARE;
simulations(12).maxRounds = 2500;
simulations(12).run(1).sources = [10 16 24 2]; % square: normal intersection, all straight through
simulations(12).run(1).targets = [6 20 4 22];
simulations(12).run(1).failed = [3 13 23];

simulations(12).run(2).sources = [16 24 10 2]; % square: normal intersection, all right turns
simulations(12).run(2).targets = [22 20 4  6];
simulations(12).run(2).failed = [3 13 23 21 25 5 1];

simulations(12).run(3).sources = [16 24 10 2]; % square: normal intersection, all left turns
simulations(12).run(3).targets = [4  6  22 20];
simulations(12).run(3).failed = [3 13 23 21 25 5 1];


simulations(13).environment = OPT_TESSELATION_SNUB_SQUARE_TILING;
simulations(13).maxRounds = 250;
simulations(13).run(1).sources = [3 4 6 5]; % snub square, 4 copies
simulations(13).run(1).targets = [10 15 24 23];% snub square, 4 copies
%targets = [10 15 11 18];% snub square, 3 copies
%sources = [3 4 6 5]; % snub square, 3 copies


simulations(14).environment = OPT_TESSELATION_PARALLELOGRAM;
simulations(14).maxRounds = 2500;
simulations(14).run(1).sources = 4; % parallelogram
simulations(14).run(1).targets = 1; % parallelogram with failures
simulations(14).run(1).failed = [3 7 11 13 9 5]; % parallelogram

simulations(15).environment = OPT_TESSELATION_PARALLELOGRAM;
simulations(15).maxRounds = 2500;
simulations(15).run(1).sources = [3 11 7 19]; % parallelogram
simulations(15).run(1).targets = [23 15 1 25]; % parallelogram with failures


simulations(16).environment = OPT_TESSELATION_SQUARE;
simulations(16).maxRounds = 2500;
simulations(16).run(1).sources = [2 24]; % one-lane bridge
simulations(16).run(1).targets = [22 4]; % one-lane bridge, squares
simulations(16).run(1).failed = [9 14 19 7 12 17   1  6 11 16 21 5 10 15 20 25]; % one-lane bridge

simulations(16).run(2).sources = [2 24]; % one-lane bridge
simulations(16).run(2).targets = [22 4]; % one-lane bridge, squares
simulations(16).run(2).failed = [9 14 19]; % one-lane bridge, dynamic failure
simulations(16).run(2).faildynamic = [7 12 17   1  6 11 16 21 5 10 15 20 25];


simulations(17).environment = OPT_TESSELATION_SQUARE;
simulations(17).maxRounds = 200;
simulations(17).run(1).sources = [5 11 21]; % source fairness assumptions
simulations(17).run(1).targets = [25 12 1];
simulations(17).run(1).failed = [];
simulations(17).run(1).faildynamic = [20 19 14 9];


simulations(18).environment = OPT_TESSELATION_SQUARE;
simulations(18).maxRounds = 500;
simulations(18).run(1).sources = [10 16]; % one-lane bridge (vertical)
simulations(18).run(1).targets = [6 20]; % one-lane bridge, (vertical)
simulations(18).run(1).failed = [1 2 3 4 5 7 8 9 17 18 19 21 22 23 24 25]; % one-lane bridge (vertical)


simulations(19).environment = OPT_TESSELATION_SQUARE;
simulations(19).maxRounds = 200;
simulations(19).run(1).sources = [2  8 11 20 24 5]; % lockcolors example
simulations(19).run(1).targets = [17 6 13 18 14 15]; % lockcolors example
simulations(19).run(1).failed = [];
simulations(19).run(1).faildynamic = [];

simulations(20).environment = OPT_TESSELATION_TRI;
simulations(20).maxRounds = 750;
simulations(20).run(1).sources = [3 11 50 40 26]; % square left-right throughput, 2 overlap (has turns)
simulations(20).run(1).targets = [23 19 25 30 41]; 
simulations(20).run(1).failed = [14 9 7 39 38 36 35];
simulations(20).run(1).faildynamic = [];

simulations(20).run(2).sources = [3 11 50 40 26]; % square left-right throughput, 2 overlap (has turns)
simulations(20).run(2).targets = [23 19 25 30 41]; 
simulations(20).run(2).failed = [];
simulations(20).run(2).faildynamic = [14 9 7 39 38 36 35];


simulations(21).environment = OPT_TESSELATION_SQUARE;
simulations(21).maxRounds = 500;
simulations(21).run(1).sources = [4 20 22 6]; % square corners
simulations(21).run(1).targets = [25 21 1 5];
simulations(21).run(1).failed = [7 8 9 12 13 14 17 18 19];


%targets = 10;  % hexagon
%sources = 16; % hexagon

%sources = [3 11]; % square left-right
%targets = [23 15]; 






%targets = [25 21 1 5 23 15]; % square corners with more
%sources = [4 20 22 6 2 16];
%failed = [7 8 9 12 13 14 17 18 19];



