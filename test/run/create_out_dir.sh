TARGET_DIR="${1:-./output}"

echo "Creating directories under: $TARGET_DIR"


datasets=(
    "citeseer"
    "dblp"
    "HPRD"
    "human"
    "maayan-figeys"
    "twitch"
    "web-Stanford"
    "wordnet-words"
    "YeastS"
    "youtube"
)

Ls=("L15" "L30" "L45" "L60")
Qs=("Q10" "Q20" "Q30" "Q40" "Q50")
Bs=("B10" "B100" "B1000" "B10000" "B100000" "B200" "B50")
# overall
for ds in "${datasets[@]}"; do
    for L in "${Ls[@]}"; do
        for Q in "${Qs[@]}"; do
            mkdir -p "$TARGET_DIR/$ds/$L/$Q"
        done
    done
    for B in "${Bs[@]}"; do
        for Q in "${Qs[@]}"; do
            mkdir -p "$TARGET_DIR/$ds/$B/$Q"
        done
    done
done

# global
global_dir="$TARGET_DIR/global/o_10_5"
mkdir -p "$global_dir"
for ds in "${datasets[@]}"; do
    for L in "${Ls[@]}"; do
        mkdir -p "$global_dir/$ds/$L"
    done
done

# degree
degrees=("D1" "D2" "D3" "D4" "D5")
Qs=("Q20" "Q30" "Q40" "Q50")
degree_dir="$TARGET_DIR/degree"
mkdir -p "$degree_dir"
for ds in "${datasets[@]}"; do
    for L in "${Ls[@]}"; do
        for D in "${degrees[@]}"; do
            for Q in "${Qs[@]}"; do
                mkdir -p "$degree_dir/$ds/$L/$D/$Q"
            done
        done
    done
done

# similar
Qs=("Q10" "Q20" "Q30" "Q40" "Q50")
Similars=("0" "25" "50" "75" "100")
similar_dir="$TARGET_DIR/similar"
mkdir -p "$similar_dir"
for Similar in "${Similars[@]}"; do
    for ds in "${datasets[@]}"; do
        for L in "${Ls[@]}"; do
            for Q in "${Qs[@]}"; do
                mkdir -p "$similar_dir/s$Similar/$ds/$L/$Q"
            done
        done
    done
done

# refine
Refines=("0" "25" "50" "75" "100")
refine_dir="$TARGET_DIR/refine"
mkdir -p "$refine_dir"
for Refine in "${Refines[@]}"; do
    for ds in "${datasets[@]}"; do
        for L in "${Ls[@]}"; do
            for Q in "${Qs[@]}"; do
                mkdir -p "$refine_dir/r$Refine/$ds/$L/$Q"
            done
        done
    done
done

# mem
Methods=("BSX" "VEQ" "GUP" "KSS" "BICE" "RM")
mem_dir="$TARGET_DIR/mem"
mkdir -p "$mem_dir"
for Method in "${Methods[@]}"; do
    for ds in "${datasets[@]}"; do
        for L in "${Ls[@]}"; do
            for Q in "${Qs[@]}"; do
                mkdir -p "$mem_dir/$Method/$ds/$L/$Q"
            done
        done
    done
done

# duplicate
datasets=(
    "citeseer"
    "Figeys"
    "YeastS"
)
Methods=("bsx" "bs1")
duplicate_dir="$TARGET_DIR/duplicate"
mkdir -p "$duplicate_dir"
for Method in "${Methods[@]}"; do
    for ds in "${datasets[@]}"; do
        mkdir -p "$duplicate_dir/$ds/$Method"
    done
done

# EvoGraph
datasets=(
    "citeseer"
    "HPRD"
    "human"
    "Figeys"
    "YeastS"
)
Scales=("5" "10" "20" "50" "100")
scale_dir="$TARGET_DIR/EvoGraph"
mkdir -p "$scale_dir"
for ds in "${datasets[@]}"; do
    for Scale in "${Scales[@]}"; do
        for Q in "${Qs[@]}"; do
            mkdir -p "$scale_dir/$ds/${ds}_${Scale}/$Q"
        done
    done
done

echo "All directories created successfully at: $TARGET_DIR"