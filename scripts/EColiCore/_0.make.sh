# # tmux session
# tmux list-sessions
# tmux attach -t 
# tmux new-session -A -s RUNNER1
# tmux new-sessioxn -A -s RUNNER2
# ps -ef | grep 'EColiCore/3.0_hit.and.down.jl'
# df -lh .

# julia -t4 --project scripts/EColiCore/2.1_gen.ecoli_core.net0.jl

# sim tools
# julia -t45 --project scripts/EColiCore/6.0_hists.jl &
# ./scripts/EColiCore/_0.make.sh

# julia -t10 --project scripts/EColiCore/2.3_gen.iJR904.net0.jl &
# for ((n=0;n<3;n++)); do
#     julia -t1 --project scripts/EColiCore/3.0_hit.and.down.jl &; 
#     julia -t1 --project scripts/EColiCore/3.1_flags.jl; 
#     julia -t1 --project scripts/EColiCore/4.0_sample.feasets.jl; 
#     julia -t1 --project scripts/EColiCore/4.1_flags.jl; 
#     julia -t1 --project scripts/EColiCore/5.0_fba.features.jl; 
#     julia -t1 --project scripts/EColiCore/5.1_flags.jl; 
# done 

# for ((n=0;n<10;n++)); do
#     julia -t1 --project scripts/EColiCore/3.0_hit.and.down.jl &
#     # julia -t1 --project scripts/EColiCore/5.0_fba.features.jl & 
# done 

for ((n=0;n<3;n++)); do
    julia -t1 --project scripts/EColiCore/3.0_hit.and.down.jl; 
    julia -t1 --project scripts/EColiCore/3.1_flags.jl; 
    julia -t1 --project scripts/EColiCore/4.0_sample.feasets.jl; 
    julia -t1 --project scripts/EColiCore/4.1_flags.jl; 
    julia -t1 --project scripts/EColiCore/5.0_fba.features.jl; 
    julia -t1 --project scripts/EColiCore/5.1_flags.jl; 
done