#source ~/.zshrc

# runs one macro and prints how long it took, compiled with ACLiC (+O) unless ACLIC= (empty), see plot.sh
run(){
    local start=$SECONDS
    root -l -b -q "${1/.cpp/.cpp${ACLIC-+O}}"
    echo ">>> $1 took $((SECONDS - start)) s"
}

mkdir -p cuts # every cut_*.root and the KE cut lists
run 'CutEvents.cpp'
run 'ReconAngleAll.cpp'
run 'CalcAPandT.cpp'
run 'CalcAPandT.cpp("fDEnergy")'
run 'CountCE.cpp'
zsh plot.sh
