set -euo pipefail

# ensure python is installed
if ! command -v python >/dev/null 2>&1
then
    echo "Please install python or activate a python environment."
    exit 1
fi

# ensure pyrosetta is installed
if ! python -c "import pyrosetta" >/dev/null 2>&1
then
    echo "Pyrosetta is not installed to active python environment; nothing to do."
    exit 1
fi

if ! python -m pybind11_stubgen >/dev/null 2>&1
then
    echo "pybind11-stubgen not installed, attempting install via pip..."
    if ! python -m pip install pybind11-stubgen
    then
        echo "pybind11-stubgen install unsuccessful, exiting. Please install pybind11-stubgen to current python environment."
    fi
fi

# run pybind11-stubgen, outputting stubs to './pyrosetta'
echo "Writing pyrosetta stubs to $(pwd)/pyrosetta..."
python -m pybind11_stubgen pyrosetta -o . || exit #make sure we make stubs correctly before copying/removing

#get site packages folder. customize if needed
sitepackages=$(python -c "import site; print(site.getsitepackages()[0])")
echo "Site-packages dir found: ${sitepackages}"

# copy stubs to site-packages
echo "Copying $(pwd)/pyrosetta/rosetta to ${sitepackages}/pyrosetta/rosetta..."
rsync -r ./pyrosetta/rosetta "${sitepackages}/pyrosetta"

#cleanup stubs
echo "Removing temporary stub directory $(pwd)/pyrosetta..."
rm -r pyrosetta

echo "Stubs copied to site-packages."
