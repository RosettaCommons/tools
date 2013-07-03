ALTER TABLE chain_endings RENAME TO chain_endings_old;
ALTER TABLE fold_trees RENAME TO fold_trees_old;
ALTER TABLE job_string_data RENAME TO job_string_data_old;
ALTER TABLE job_string_real_data RENAME TO job_string_real_data_old;
ALTER TABLE job_string_string_data RENAME TO job_string_string_data_old;
ALTER TABLE jumps RENAME TO jumps_old;
ALTER TABLE nonprotein_residue_angles RENAME TO nonprotein_residue_angles_old;
ALTER TABLE nonprotein_residue_conformation RENAME TO nonprotein_residue_conformation_old;
ALTER TABLE pose_comments RENAME TO pose_comments_old;
ALTER TABLE pose_conformations RENAME TO pose_conformations_old;
ALTER TABLE protein_residue_conformation RENAME TO protein_residue_conformation_old;
ALTER TABLE protocols RENAME TO protocols_old;
ALTER TABLE residue_atom_coords RENAME TO residue_atom_coords_old;
ALTER TABLE residue_pdb_confidence RENAME TO residue_pdb_confidence_old;
ALTER TABLE residue_pdb_identification RENAME TO residue_pdb_identification_old;
ALTER TABLE residues RENAME TO residues_old;
ALTER TABLE structure_scores RENAME TO structure_scores_old;
ALTER TABLE structures RENAME TO structures_old;


CREATE TABLE struct_id_uuid_lookup (
  struct_id_int INTEGER PRIMARY KEY AUTOINCREMENT,
  struct_id_uuid BLOB NOT NULL);

INSERT INTO struct_id_uuid_lookup (struct_id_uuid) SELECT struct_id FROM structures_old;

--CREATE VIEW name_view AS
--  SELECT 
--    struct_id_int as struct_id,
--    remaining_column
--  FROM
--    name_old JOIN struct_id_uuid_lookup ON name_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;
    

CREATE VIEW chain_endings_view AS
  SELECT
    struct_id_int as struct_id,
    end_pos
  FROM
    chain_endings_old JOIN struct_id_uuid_lookup ON chain_endings_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;
CREATE VIEW fold_trees_view AS
  SELECT
    struct_id_int as struct_id,
    start_res,
    start_atom,
    stop_res,
    stop_atom,
    label,
    keep_stub_in_residue
  FROM
    fold_trees_old JOIN struct_id_uuid_lookup ON fold_trees_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW job_string_data_view AS
  SELECT
    struct_id_int as struct_id,
    data_key
  FROM
    job_string_data_old JOIN struct_id_uuid_lookup ON job_string_data_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW job_string_real_data_view AS
  SELECT
    struct_id_int as struct_id,
    data_key,
    data_value
  FROM
    job_string_real_data_old JOIN struct_id_uuid_lookup ON job_string_real_data_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW job_string_string_data_view AS
  SELECT
    struct_id_int as struct_id,
    data_key,
    data_value
  FROM
    job_string_string_data_old JOIN struct_id_uuid_lookup ON job_string_string_data_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW jumps_view AS
  SELECT
    struct_id_int as struct_id,
    jump_id,
    xx,
    xy,
    xz,
    yx,
    yy,
    yz,
    zx,
    zy,
    zz,
    x,
    y,
    z
  FROM
    jumps_old JOIN struct_id_uuid_lookup ON jumps_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW nonprotein_residue_angles_view AS
  SELECT
    struct_id_int as struct_id,
    seqpos,
    chinum,
    chiangle
  FROM
    nonprotein_residue_angles_old JOIN struct_id_uuid_lookup ON nonprotein_residue_angles_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW nonprotein_residue_conformation_view AS
  SELECT
    struct_id_int as struct_id,
    seqpos,
    phi,
    psi,
    omega
  FROM
    nonprotein_residue_conformation_old JOIN struct_id_uuid_lookup ON nonprotein_residue_conformation_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW pose_comments_view AS
  SELECT
    struct_id_int as struct_id,
    comment_key,
    value
  FROM
    pose_comments_old JOIN struct_id_uuid_lookup ON pose_comments_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW pose_conformations_view AS
  SELECT
    struct_id_int as struct_id,
    annotated_sequence,
    total_residue,
    fullatom
  FROM
    pose_conformations_old JOIN struct_id_uuid_lookup ON pose_conformations_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW protein_residue_conformation_view AS
  SELECT
    struct_id_int as struct_id,
    seqpos,
    secstruct,
    phi,
    psi,
    omega,
    chi1,
    chi2,
    chi3,
    chi4
  FROM
    protein_residue_conformation_old JOIN struct_id_uuid_lookup ON protein_residue_conformation_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW residue_atom_coords_view AS
  SELECT
    struct_id_int as struct_id,
    seqpos,
    atomno,
    x,
    y,
    z
  FROM
    residue_atom_coords_old JOIN struct_id_uuid_lookup ON residue_atom_coords_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW residue_pdb_confidence_view AS
  SELECT
    struct_id_int as struct_id,
    residue_number,
    max_temperature,
    max_bb_temperature,
    max_sc_temperature,
    min_occupancy,
    min_bb_occupancy,
    min_sc_occupancy
  FROM
    residue_pdb_confidence_old JOIN struct_id_uuid_lookup ON residue_pdb_confidence_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW residue_pdb_identification_view AS
  SELECT
    struct_id_int as struct_id,
    residue_number,
    chain_id,
    insertion_code,
    pdb_residue_number
  FROM
    residue_pdb_identification_old JOIN struct_id_uuid_lookup ON residue_pdb_identification_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW residues_view AS
  SELECT
    struct_id_int as struct_id,
    resNum,
    name3,
    res_type
  FROM
    residues_old JOIN struct_id_uuid_lookup ON residues_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW structure_scores_view AS
  SELECT
    batch_id,
    struct_id_int as struct_id,
    score_type_id,
    score_value
  FROM
    structure_scores_old JOIN struct_id_uuid_lookup ON structure_scores_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;

CREATE VIEW structures_view AS
  SELECT
    struct_id_int as struct_id,
    batch_id,
    tag,
    input_tag
FROM
  structures_old JOIN struct_id_uuid_lookup ON structures_old.struct_id = struct_id_uuid_lookup.struct_id_uuid;
