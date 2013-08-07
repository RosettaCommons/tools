CREATE TABLE chain_endings(
	struct_id INTEGER NOT NULL,
	end_pos INTEGER,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED);

CREATE TABLE fold_trees(
	struct_id INTEGER NOT NULL,
	start_res INTEGER,
	start_atom TEXT,
	stop_res INTEGER,
	stop_atom TEXT,
	label INTEGER,
	keep_stub_in_residue INTEGER,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED);

CREATE TABLE job_string_data(
	struct_id INTEGER NOT NULL,
	data_key TEXT,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id, data_key));

CREATE TABLE job_string_real_data(
	struct_id INTEGER NOT NULL,
	data_key TEXT,
	data_value REAL,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id, data_key));

CREATE TABLE job_string_string_data(
	struct_id INTEGER NOT NULL,
	data_key TEXT,
	data_value TEXT,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id, data_key));

CREATE TABLE jumps(
	struct_id INTEGER NOT NULL,
	jump_id INTEGER,
	xx REAL,
	xy REAL,
	xz REAL,
	yx REAL,
	yy REAL,
	yz REAL,
	zx REAL,
	zy REAL,
	zz REAL,
	x REAL,
	y REAL,
	z REAL,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED);
CREATE TABLE nonprotein_residue_angles(
	struct_id INTEGER NOT NULL,
	seqpos INTEGER NOT NULL,
	chinum INTEGER NOT NULL,
	chiangle REAL NOT NULL,
	FOREIGN KEY (struct_id, seqpos) REFERENCES residues(struct_id, resNum) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id, seqpos, chinum));
CREATE TABLE nonprotein_residue_conformation(
	struct_id INTEGER NOT NULL,
	seqpos INTEGER NOT NULL,
	phi REAL NOT NULL,
	psi REAL NOT NULL,
	omega REAL NOT NULL,
	FOREIGN KEY (struct_id, seqpos) REFERENCES residues(struct_id, resNum) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id, seqpos));
CREATE TABLE pose_comments(
	struct_id INTEGER NOT NULL,
	comment_key TEXT NOT NULL,
	value TEXT NOT NULL,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id, comment_key));
CREATE TABLE pose_conformations(
	struct_id INTEGER NOT NULL,
	annotated_sequence TEXT,
	total_residue INTEGER,
	fullatom INTEGER,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED);
CREATE TABLE protein_residue_conformation(
	struct_id INTEGER NOT NULL,
	seqpos INTEGER NOT NULL,
	secstruct TEXT NOT NULL,
	phi REAL NOT NULL,
	psi REAL NOT NULL,
	omega REAL NOT NULL,
	chi1 REAL NOT NULL,
	chi2 REAL NOT NULL,
	chi3 REAL NOT NULL,
	chi4 REAL NOT NULL,
	FOREIGN KEY (struct_id, seqpos) REFERENCES residues(struct_id, resNum) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id, seqpos));
CREATE TABLE residue_atom_coords(
	struct_id INTEGER NOT NULL,
	seqpos INTEGER NOT NULL,
	atomno INTEGER NOT NULL,
	x REAL NOT NULL,
	y REAL NOT NULL,
	z REAL NOT NULL,
	FOREIGN KEY (struct_id, seqpos) REFERENCES residues(struct_id, resNum) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id, seqpos, atomno));
CREATE TABLE residue_pdb_confidence(
	struct_id INTEGER NOT NULL,
	residue_number INTEGER NOT NULL,
	max_temperature REAL NOT NULL,
	max_bb_temperature REAL NOT NULL,
	max_sc_temperature REAL NOT NULL,
	min_occupancy REAL NOT NULL,
	min_bb_occupancy REAL NOT NULL,
	min_sc_occupancy REAL NOT NULL,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id, residue_number));
CREATE TABLE residue_pdb_identification(
	struct_id INTEGER NOT NULL,
	residue_number INTEGER NOT NULL,
	chain_id TEXT NOT NULL,
	insertion_code TEXT NOT NULL,
	pdb_residue_number INTEGER NOT NULL,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id, residue_number));
CREATE TABLE residues(
	struct_id INTEGER NOT NULL,
	resNum INTEGER NOT NULL,
	name3 TEXT NOT NULL,
	res_type TEXT NOT NULL,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (struct_id, resNum));
CREATE TABLE structure_scores(
	batch_id INTEGER NOT NULL,
	struct_id INTEGER NOT NULL,
	score_type_id INTEGER NOT NULL,
	score_value INTEGER NOT NULL,
	FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY (batch_id, score_type_id) REFERENCES score_types(batch_id, score_type_id) DEFERRABLE INITIALLY DEFERRED,
	PRIMARY KEY (batch_id, struct_id, score_type_id));
CREATE TABLE structures(
	struct_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
	batch_id INTEGER,
	tag TEXT,
	input_tag TEXT,
	FOREIGN KEY (batch_id) REFERENCES batches(batch_id) DEFERRABLE INITIALLY DEFERRED);

INSERT INTO chain_endings SELECT * FROM chain_endings_view;
INSERT INTO fold_trees SELECT * FROM fold_trees_view;
INSERT INTO job_string_data SELECT * FROM job_string_data_view;
INSERT INTO job_string_real_data SELECT * FROM job_string_real_data_view;
INSERT INTO job_string_string_data SELECT * FROM job_string_string_data_view;
INSERT INTO jumps SELECT * FROM jumps_view;
INSERT INTO nonprotein_residue_angles SELECT * FROM nonprotein_residue_angles_view;
INSERT INTO nonprotein_residue_conformation SELECT * FROM nonprotein_residue_conformation_view;
INSERT INTO pose_comments SELECT * FROM pose_comments_view;
INSERT INTO pose_conformations SELECT * FROM pose_conformations_view;
INSERT INTO protein_residue_conformation SELECT * FROM protein_residue_conformation_view;
INSERT INTO residue_atom_coords SELECT * FROM residue_atom_coords_view;
INSERT INTO residue_pdb_confidence SELECT * FROM residue_pdb_confidence_view;
INSERT INTO residue_pdb_identification SELECT * FROM residue_pdb_identification_view;
INSERT INTO residues SELECT * FROM residues_view;
INSERT INTO structure_scores SELECT * FROM structure_scores_view;
INSERT INTO structures SELECT * FROM structures_view;
