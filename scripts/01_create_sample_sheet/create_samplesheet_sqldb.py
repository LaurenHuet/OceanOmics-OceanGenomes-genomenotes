#!/usr/bin/env python3
# singularity run $SING/psycopg2:0.1.sif python create_samplesheet_sqldb.py ~/postgresql_details/oceanomics.cfg
import sys
import os
import getpass
from pathlib import Path
import configparser
from datetime import date
import pandas as pd
import psycopg2
from psycopg2.extras import execute_values

def load_db_config(config_file: str) -> dict:
    p = Path(config_file)
    if not p.exists():
        raise FileNotFoundError(f"❌ Config file '{config_file}' does not exist.")
    cfg = configparser.ConfigParser()
    cfg.read(config_file)
    if 'postgres' not in cfg:
        raise ValueError("❌ Missing [postgres] section in config file.")
    required = ['dbname', 'user', 'password', 'host', 'port']
    missing = [k for k in required if not cfg.has_option('postgres', k)]
    if missing:
        raise ValueError(f"❌ Missing keys in [postgres]: {missing}")
    return {
        'dbname': cfg.get('postgres', 'dbname'),
        'user': cfg.get('postgres', 'user'),
        'password': cfg.get('postgres', 'password'),
        'host': cfg.get('postgres', 'host'),
        'port': cfg.getint('postgres', 'port')
    }

# =====================================
# OG IDs for the samplesheet (modify as needed)
# =====================================
og_ids = [
    'OG38'
]

# =====================================
# SQL function definition template
# - assembly path
# - busco_genes path
# - bioproject_id pulled from sample.ncbi_bioproject_id_lvl_3_hifi
# - tolid pulled from sample.tol_id
# =====================================
create_function_sql_template = """
CREATE OR REPLACE FUNCTION build_nfcore_samplesheet_rows(in_og_ids text[])
RETURNS TABLE (
  sample        text,
  hifi_dir      text,
  hic_dir       text,
  assembly      text,
  busco_genes   text,
  bioproject_id text,
  version       text,
  date          text,
  tolid         text,
  taxid         bigint,
  species       text
)
LANGUAGE sql
AS $$
WITH p AS (
  SELECT unnest(in_og_ids) AS og_id
),
latest_seq AS (
  SELECT DISTINCT ON (seq.og_id)
         seq.og_id,
         seq.seq_date::date AS seq_date
  FROM sequencing seq
  JOIN p ON seq.og_id = p.og_id
  WHERE seq.technology = 'PacBio'
  ORDER BY seq.og_id, seq.seq_date DESC
),
smp AS (
  SELECT DISTINCT ON (s.og_id)
         s.og_id,
         s.nominal_species_id,
         s.tol_id,
         s.ncbi_bioproject_id_lvl_3_hifi
  FROM sample s
  JOIN p ON s.og_id = p.og_id
  ORDER BY s.og_id
)
SELECT DISTINCT ON (p.og_id)
  p.og_id AS sample,
  '/scratch/pawsey0964/$USER/genomenotes/'||p.og_id||'/hifi' AS hifi_dir,
  '/scratch/pawsey0964/$USER/genomenotes/'||p.og_id||'/hic'  AS hic_dir,
  '/scratch/pawsey0964/$USER/genomenotes/'||p.og_id||'/assembly' AS assembly,
  '/scratch/pawsey0964/$USER/genomenotes/'||p.og_id||'/busco' AS busco_genes,
  smp.ncbi_bioproject_id_lvl_3_hifi AS bioproject_id,
  CASE WHEN rg.og_id IS NOT NULL THEN 'hic2' ELSE 'hic1' END AS version,
  CASE WHEN ls.seq_date IS NOT NULL THEN 'v'||to_char(ls.seq_date,'YYMMDD') END AS date,
  smp.tol_id AS tolid,
  sp.ncbi_taxon_id AS taxid,
  sp.species
FROM p
LEFT JOIN ref_genomes rg ON rg.og_id = p.og_id
LEFT JOIN latest_seq ls  ON ls.og_id = p.og_id
LEFT JOIN smp ON smp.og_id = p.og_id
LEFT JOIN species sp ON sp.species = smp.nominal_species_id
ORDER BY p.og_id;
$$;
"""

# get runtime user in a portable way and sanitize
runner_user = os.environ.get('USER') or getpass.getuser()
runner_user = runner_user.replace("'", "").replace("/", "")

# replace $USER placeholder in the template
create_function_sql = create_function_sql_template.replace("$USER", runner_user)

# =====================================
# Main: connect, create function, call it, save CSV
# =====================================
def main():
    if len(sys.argv) < 2:
        print("Usage: create_samplesheet_sqldb.py <db_config_file>")
        sys.exit(1)

    config_file = sys.argv[1]
    conn = None
    cur = None
    try:
        db_params = load_db_config(config_file)
        conn = psycopg2.connect(**db_params)
        cur = conn.cursor()

        # DROP existing function first (can't change OUT columns with CREATE OR REPLACE)
        cur.execute("DROP FUNCTION IF EXISTS build_nfcore_samplesheet_rows(text[]);")
        conn.commit()

        # Create the SQL function (with busco_genes and bioproject_id)
        cur.execute(create_function_sql)
        conn.commit()

        # Call the function with OG list
        query = "SELECT * FROM build_nfcore_samplesheet_rows(%s);"
        cur.execute(query, (og_ids,))

        rows = cur.fetchall()
        cols = [desc[0] for desc in cur.description]
        df = pd.DataFrame(rows, columns=cols)

        # Ensure taxid is nullable integer type (remove trailing .0 if present)
        if 'taxid' in df.columns:
            try:
                df['taxid'] = df['taxid'].astype('Int64')
            except Exception:
                df['taxid'] = pd.to_numeric(df['taxid'], errors='coerce').astype('Int64')

        # Identify and print rows with missing values
        missing_rows = df[df.isnull().any(axis=1)]
        if not missing_rows.empty:
            print("\nRows with missing values:\n")
            print(missing_rows)

        # Save CSV to current working directory
        today_str = date.today().strftime("%Y%m%d")
        current_dir = os.getcwd()
        output_path = os.path.join(current_dir, f"samplesheet_{today_str}.csv")
        df.to_csv(output_path, index=False)

        print(f"Samplesheet saved to: {output_path}")

    finally:
        if cur is not None:
            cur.close()
        if conn is not None:
            conn.close()

if __name__ == "__main__":
    main()