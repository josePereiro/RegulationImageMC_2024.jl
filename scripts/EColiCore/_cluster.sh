# tmux session
tmux list-sessions
tmux attach -t 
tmux new-session -A -s RUNNER1
tmux new-session -A -s RUNNER2


# sim tools
julia -t10 --project scripts/EColiCore/2.3_gen.iJR904.net0.jl &
julia -t1 --project scripts/EColiCore/3.0_hit.and.down.jl &
julia -t1 --project scripts/EColiCore/3.1_flag.duplicates &
julia -t1 --project scripts/EColiCore/4.0_sample.feasets.jl &
julia -t1 --project scripts/EColiCore/4.1_flag.duplicates &
julia -t1 --project scripts/EColiCore/5.0_fba.features.jl &
julia -t1 --project scripts/EColiCore/6.0_hists.jl &