for i in "$@"; do 
    sed -i  's/Bin_Screener/BinScreener/g' $i
    sed -i  's/Base_Sugar_Rotamer/BaseSugarRotamer/g' $i
    sed -i  's/Base_Sampler_Util/BaseSamplerUtil/g' $i
    sed -i  's/Parameters_Setup/ParametersSetup/g' $i
    sed -i  's/Generator_Wrapper/GeneratorWrapper/g' $i
done