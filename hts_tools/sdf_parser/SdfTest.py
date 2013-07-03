import SdfFile

if __name__ == "__main__":
    sdf_file = SdfFile.SdfFile()
    sdf_file.add_file("test_file.sdf")
    test_data = ["a","b"]
    sdf_file.get_record("ZINC10504980").add_data_entry("TestData",test_data)
    sdf_file.write_sdf_file("out_file.sdf")
    
    print "testing generator"
    for molfile in SdfFile.parse_sdf("test_file.sdf"):
        print molfile.get_molecule_name()