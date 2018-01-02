# test small chain

# model = Tmp.chain(4, :periodic, 0., complex(1.0))
# Tmp.runVMC(model, 2, true)

# test VMC for chain

# model = Tmp.chain(32, :periodic, 0., complex(1.0))
# Tmp.runVMC(model, 10000)

# test VMC for kagome strip

# list = [
#     0.5 -1.2; 0.5 -1.4; 0.5 -1.7; 0.5 -2.1;
#     0.6 -1.3; 0.6 -1.5; 0.6 -1.9; 0.6 -2.5;
#     0.7 -1.4; 0.7 -1.7; 0.7 -2.1; 0.7 -3.1;
#     0.8 -1.5; 0.8 -2.0; 0.8 -2.3; 0.8 -3.6;
#     0.9 -1.6; 0.9 -2.1; 0.9 -2.7; 0.9 -4.1;
#     1.0 -1.8; 1.0 -2.4; 1.0 -3.1; 1.0 -4.8;
#     1.1 -1.9; 1.1 -2.6; 1.1 -3.4; 1.1 -5.6;
#     1.2 -2.1; 1.2 -2.9; 1.2 -3.9; 1.2 -6.4;
#     1.3 -2.3; 1.3 -3.2; 1.3 -4.3; 1.3 -7.0;
# ]

list = [
    0.9 -1.6; 0.9 -2.1; 0.9 -2.7; 0.9 -4.1;
    1.0 -1.8; 1.0 -2.4; 1.0 -3.1; 1.0 -4.8;
]

# list = [
#     1.0 -2.5;
# ]

for i=1:size(list)[1]
    t = list[i, 1]
    m = list[i, 2]
    model = Tmp.kagomestrip_LC(64, :open, m, complex(t))
    Tmp.runVMC(model, 16000000)
end
