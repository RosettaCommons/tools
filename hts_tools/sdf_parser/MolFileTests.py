import MolFile
import unittest

        
class TestParserFunctions(unittest.TestCase):
    
    def test_convert_or_blank(self):
        self.assertEqual(MolFile.convert_or_blank("123",int),123)
        self.assertEqual(MolFile.convert_or_blank("   ",int),"")
    
    def test_parse_counts(self):
        count_line = "111222333444555666777888999111999 V2000"
        count_data = MolFile.parse_counts(count_line)
        self.assertEqual(count_data.n_atoms,111)
        self.assertEqual(count_data.n_bonds,222)
        self.assertEqual(count_data.n_lists,333)
        self.assertEqual(count_data.unused_a,444)
        self.assertEqual(count_data.chiral_flag,555)
        self.assertEqual(count_data.n_stext,666)
        self.assertEqual(count_data.unused_b,777)
        self.assertEqual(count_data.unused_c,888)
        self.assertEqual(count_data.unused_d,999)
        self.assertEqual(count_data.unused_e,111)
        self.assertEqual(count_data.n_prop,999)
        self.assertEqual(count_data.version," V2000")

    def test_output_counts(self):
        count_line = "111222333   555666777888999111999 V2000"
        count_data = MolFile.parse_counts(count_line)
        new_count_line = MolFile.counts_to_string(count_data)
        self.assertEqual(new_count_line,count_line)

    def test_parse_atom(self):
        atom_line = "11115.5159-2223.274233333.7439 CCC11222333444555666777888999110111112"
        atom_data = MolFile.parse_atom(atom_line,1)
        self.assertEqual(atom_data.x_coord,11115.5159)
        self.assertEqual(atom_data.y_coord,-2223.2742)
        self.assertEqual(atom_data.z_coord,33333.7439)
        self.assertEqual(atom_data.symbol,"CCC")
        self.assertEqual(atom_data.mass_diff,11)
        self.assertEqual(atom_data.charge,222)
        self.assertEqual(atom_data.stereo_parity,333)
        self.assertEqual(atom_data.hydrogen_count,444)
        self.assertEqual(atom_data.stereo_care,555)
        self.assertEqual(atom_data.valence,666)
        self.assertEqual(atom_data.h0_designator,777)
        self.assertEqual(atom_data.unused_a,888)
        self.assertEqual(atom_data.unused_b,999)
        self.assertEqual(atom_data.atom_atom_map,110)
        self.assertEqual(atom_data.invert_flag,111)
        self.assertEqual(atom_data.exact_change,112)
        self.assertEqual(atom_data.index,1)

    def test_output_atom(self):
        atom_line = "11115.5159-2223.274233333.7439 CCC11222333444555666777888999110111112"
        atom_data = MolFile.parse_atom(atom_line,1)
        new_atom_line = MolFile.atom_to_string(atom_data)
        self.assertEqual(new_atom_line,atom_line)

    def test_parse_bond(self):
        bond_line = "111222333444555666777"
        bond_data = MolFile.parse_bond(bond_line)
        self.assertEqual(bond_data.first_atom,111)
        self.assertEqual(bond_data.second_atom,222)
        self.assertEqual(bond_data.bond_type,333)
        self.assertEqual(bond_data.bond_stereo,444)
        self.assertEqual(bond_data.unused_a,555)
        self.assertEqual(bond_data.bond_topology,666)
        self.assertEqual(bond_data.reaction_center,777)

    def test_output_bond(self):
        bond_line = "111222333444555666777"
        bond_data = MolFile.parse_bond(bond_line)
        new_bond_line = MolFile.bond_to_string(bond_data)
        self.assertEqual(new_bond_line,bond_line)
        
    def test_parse_header(self):
        
        header_lines = [
        "ZINC00000221",
        "  -OEChem-08230709433D",
        "comments"
        ]
        
        header_data = MolFile.parse_header(header_lines)
        self.assertEqual(header_data.name,"ZINC00000221")
        self.assertEqual(header_data.program,"  -OEChem-08230709433D")
        self.assertEqual(header_data.comments,"comments")
        
    def test_output_header(self):
        header_lines = [
        "ZINC00000221",
        "  -OEChem-08230709433D",
        "comments"
        ]
        header_data = MolFile.parse_header(header_lines)
        new_header_lines = MolFile.header_to_strings(header_data)
        self.assertEqual(new_header_lines,header_lines)
        
    def test_parse_data(self):
        data_lines = [
            ">  <ChemCompId>",
            "FVF",
            "",
            ">  <PdbId>",
            "7CPA",
            "",
            ">  <ChainId>",
            "A",
            "",
            ">  <ResidueNumber>",
            "309",
            "",
            ">  <InsertionCode>",
            "",
            ""
            ]
        
        data_block = MolFile.parse_data_block(data_lines)
        self.assertEqual(data_block["ChemCompId"],["FVF"])
        self.assertEqual(data_block["PdbId"],["7CPA"])
        self.assertEqual(data_block["ChainId"],["A"])
        self.assertEqual(data_block["ResidueNumber"],["309"])
        self.assertEqual(data_block["InsertionCode"],[""])
        
    def test_output_data(self):
        data_lines = [
            ">  <ChemCompId>",
            "FVF",
            "",
            ">  <PdbId>",
            "7CPA",
            "",
            ">  <ChainId>",
            "A",
            "",
            ">  <ResidueNumber>",
            "309",
            "",
            ">  <InsertionCode>",
            "",
            ""
            ]
        data_block = MolFile.parse_data_block(data_lines)
        data_strings = MolFile.data_block_to_strings(data_block)
        correct_answer = [
            '>  <PdbId>', 
            '7CPA', 
            '', 
            '>  <InsertionCode>', 
            '', 
            '', 
            '>  <ChainId>', 
            'A', 
            '', 
            '>  <ResidueNumber>', 
            '309', 
            '', 
            '>  <ChemCompId>', 
            'FVF', 
            '']
        self.assertEqual(correct_answer,data_strings)
        
        
class TestMolFile(unittest.TestCase):
    def test_make_molfile(self):
        molfile_lines = [
            "ZINC00000221",
            "  -OEChem-08230709433D",
            " line 3",
            " 33 34  0     1  0  0  0  0  0999 V2000",
            "   -5.5159   -3.2742    3.7439 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -5.6927   -3.4595    2.3015 N   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -6.7971   -3.9197    1.6820 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -5.3635   -3.5216    0.1113 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -4.7590   -3.1968    1.3192 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -3.4334   -2.7116    1.3164 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -2.8803   -2.4244    2.3628 O   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -2.7877   -2.5661    0.1406 N   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -3.3846   -2.8846   -1.0208 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -2.7701   -2.7448   -2.0608 O   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -4.6458   -3.3542   -1.0619 N   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -5.2599   -3.6892   -2.3491 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -1.4133   -2.0590    0.1302 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -1.4333   -0.5336    0.0129 C   0  0  1  0  0  0  0  0  0  0  0  0",
            "   -1.9701   -0.1109    0.8621 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "    0.0021   -0.0041    0.0020 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -0.0173    1.4248    0.0099 O   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -2.0889   -0.1578   -1.2000 O   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -5.1169   -4.1892    4.1816 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -4.8215   -2.4532    3.9223 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -6.4779   -3.0422    4.2011 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -7.7099   -4.2162    2.1773 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -4.5494   -3.4940   -3.1524 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -5.5360   -4.7436   -2.3565 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -6.1512   -3.0792   -2.4962 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -0.8764   -2.4817   -0.7190 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -0.9132   -2.3457    1.0553 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "    0.5123   -0.3556   -0.8948 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "    0.5293   -0.3651    0.8851 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "    0.8606    1.8301    0.0037 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -1.6646   -0.4999   -1.9987 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "   -6.5975   -3.9546    0.3879 N   0  3  0  0  0  0  0  0  0  0  0  0",
            "   -7.2846   -4.2669   -0.3083 H   0  0  0  0  0  0  0  0  0  0  0  0",
            "  1  2  1  0  0  0  0",
            "  1 19  1  0  0  0  0",
            "  1 20  1  0  0  0  0",
            "  1 21  1  0  0  0  0",
            "  2  5  1  0  0  0  0",
            "  2  3  1  0  0  0  0",
            "  3 22  1  0  0  0  0",
            "  3 32  2  0  0  0  0",
            "  4 11  1  0  0  0  0",
            "  4  5  2  0  0  0  0",
            "  4 32  1  0  0  0  0",
            "  5  6  1  0  0  0  0",
            "  6  7  2  0  0  0  0",
            "  6  8  1  0  0  0  0",
            "  8  9  1  0  0  0  0",
            "  8 13  1  0  0  0  0",
            "  9 10  2  0  0  0  0",
            "  9 11  1  0  0  0  0",
            " 11 12  1  0  0  0  0",
            " 12 23  1  0  0  0  0",
            " 12 24  1  0  0  0  0",
            " 12 25  1  0  0  0  0",
            " 13 14  1  0  0  0  0",
            " 13 26  1  0  0  0  0",
            " 13 27  1  0  0  0  0",
            " 14 15  1  0  0  0  0",
            " 14 16  1  0  0  0  0",
            " 14 18  1  0  0  0  0",
            " 16 17  1  0  0  0  0",
            " 16 28  1  0  0  0  0",
            " 16 29  1  0  0  0  0",
            " 17 30  1  0  0  0  0",
            " 18 31  1  0  0  0  0",
            " 32 33  1  0  0  0  0",
            "M  CHG  1  32   1 ",
            "M  END",
            ">  <ChainId>",
            "A",
            "",
            ">  <ResidueNumber>",
            "309",
            ""
        ]
        mol_test = MolFile.MolFile(molfile_lines)
        self.assertEqual(mol_test.data_to_string_list(),molfile_lines)
        
        test_data = ["a","b"]
        mol_test.add_data_entry("TestEntry",test_data)
        correct_answer = [
        '>  <TestEntry>', 
        'a', 
        'b', 
        '', 
        '>  <ChainId>', 
        'A', 
        '', 
        '>  <ResidueNumber>', 
        '309', 
        '']

        self.assertEqual(MolFile.data_block_to_strings(mol_test.data_map),correct_answer)
        
if __name__ == '__main__':
    unittest.main()