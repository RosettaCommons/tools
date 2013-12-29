find ./ -name '*.??' | xargs sed -i '~' 's/protocols\/rna/protocols\/farna/g'
find ./ -name '*.??' | xargs sed -i '~' 's/protocols::rna/protocols::farna/g'
cleanup
