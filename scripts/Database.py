#!/usr/bin/python
import sqlite3
import os
import re
import traceback


from . import myUtil
from . import ParseReports
from . import Csb_finder

"""
This module keeps all routines for database input and working routines on the database
if looking for output routines invoked by the user got to Output.py module
"""

########## Write output to Database Routines ##########

def create_database(database):
    """
    11.9.22
    Creating a database file for midterm storage of results
    
    Args:
        database    Pathway to database file
    """
    #create database
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute('''CREATE TABLE Genomes (
        genomeID    varchar(32)     PRIMARY KEY     NOT NULL,
        Superkingdom      varchar(128)              DEFAULT 'NULL',
        Clade       varchar(128)                    DEFAULT 'NULL',
        Phylum      varchar(128)                    DEFAULT 'NULL',
        Class       varchar(128)                    DEFAULT 'NULL',
        Ordnung     varchar(128)                    DEFAULT 'NULL',
        Family      varchar(128)                    DEFAULT 'NULL',
        Genus       varchar(128)                    DEFAULT 'NULL',
        Species     varchar(128)                    DEFAULT 'NULL',
        Strain      varchar(128)                    DEFAULT 'NULL',
        TypeStrain  tinyint(4)                      DEFAULT NULL,
        Completeness decimal(5,2)                   DEFAULT NULL,
        Contamination decimal(5,2)                  DEFAULT NULL,
        dRep        tinyint(1)                      DEFAULT NULL,
        NCBITaxon   int(11)                         DEFAULT NULL,
        NCBIProject int(11)                         DEFAULT NULL,
        NCBIBioproject varchar(32)                  DEFAULT NULL,
        NCBIBiosample varchar(32)                   DEFAULT NULL,
        NCBIAssembly varchar(32)                    DEFAULT NULL
        );''')
        
        cur.execute('''CREATE TABLE Clusters (
        clusterID   varchar(32)     PRIMARY KEY     NOT NULL,
        genomeID    varchar(32)                     NOT NULL,
        CONSTRAINT fk_genomeID FOREIGN KEY (genomeID) REFERENCES Genomes(genomeID) ON DELETE CASCADE ON UPDATE CASCADE
        );''')
        
        cur.execute('''CREATE TABLE Keywords (
        ID          integer         PRIMARY KEY     AUTOINCREMENT,
        clusterID   varchar(32)                     NOT NULL,
        keyword     varchar(32)                     NOT NULL,
        completeness    varchar(32)                 DEFAULT NULL,
        collinearity    varchar(32)                 DEFAULT NULL,
        CONSTRAINT fk_clusterID FOREIGN KEY (clusterID) REFERENCES Clusters(clusterID) ON DELETE CASCADE ON UPDATE CASCADE
        );''')
        
        cur.execute('''CREATE TABLE Proteins (
        proteinID   varchar(128)     PRIMARY KEY     NOT NULL,
        genomeID    varchar(64)                     NOT NULL,
        clusterID   varchar(64)                     DEFAULT NULL, 
        locustag    varchar(32)                     DEFAULT NULL,
        contig      varchar(32)                     DEFAULT NULL,
        start       int(11)                         DEFAULT NULL,
        end         int(11)                         DEFAULT NULL,
        strand      varchar(1)                      DEFAULT NULL,
        dom_count   smallint(6)                     DEFAULT NULL,
        sequence    varchar(4096)                   DEFAULT NULL,
        UNIQUE(proteinID,genomeID),
        CONSTRAINT fk_genomeID FOREIGN KEY (genomeID) REFERENCES Genomes(genomeID) ON DELETE CASCADE ON UPDATE CASCADE,
        CONSTRAINT fk_clusterID FOREIGN KEY (clusterID) REFERENCES Clusters(clusterID) ON DELETE SET NULL ON UPDATE CASCADE
        );''')
        
        cur.execute('''CREATE TABLE Domains (
        ID          integer         PRIMARY KEY     AUTOINCREMENT,
        proteinID   varchar(32)                     NOT NULL,
        domain      varchar(32)                     DEFAULT NULL,
        score       smallint(6)                     DEFAULT NULL,
        domStart    int(11)                             DEFAULT NULL,
        domEnd      int(11)                             DEFAULT NULL,
        CONSTRAINT fk_proteinID FOREIGN KEY (proteinID) REFERENCES Proteins(proteinID) ON DELETE CASCADE ON UPDATE CASCADE
        );''')
        
              
    con.commit()
    con.close()
        #res = cur.execute("Select name from sqlite_master")
        #print(res.fetchall())
    
        
    return 
    
def index_database(database):
    """
    12.11.22
    Indexes the relevant columns in the database to improve search and join performance.
    This is especially important for large databases (> 3000 genomes).
    """
    print(f"Indexing database {database}", end="\r")
    
    # List of indexes to create, each represented as a tuple (index_name, create_statement)
    indexes = [
        ("tab_dom_pid_index", "CREATE INDEX IF NOT EXISTS tab_dom_pid_index ON Domains(proteinID)"),
        ("tab_prot_gid_index", "CREATE INDEX IF NOT EXISTS tab_prot_gid_index ON Proteins(genomeID)"),
        ("tab_prot_cid_index", "CREATE INDEX IF NOT EXISTS tab_prot_cid_index ON Proteins(clusterID)"),
        ("tab_key_cid_index", "CREATE INDEX IF NOT EXISTS tab_key_cid_index ON Keywords(clusterID)"),
        ("tab_clus_gid_index", "CREATE INDEX IF NOT EXISTS tab_clus_gid_index ON Clusters(genomeID)"),
        ("tab_dom_did_index", "CREATE INDEX IF NOT EXISTS tab_dom_did_index ON Domains(domain)"),
        ("tab_key_kid_index", "CREATE INDEX IF NOT EXISTS tab_key_kid_index ON Keywords(keyword)")
    ]
    
    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            cur.execute("""PRAGMA foreign_keys = ON;""")
            
            # Start a transaction for batch index creation
            with con:
                for index_name, create_stmt in indexes:
                    cur.execute(create_stmt)
                
    except sqlite3.Error as e:
        print(f"An error occurred while indexing the database: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    print(f"Indexing database {database} -- ok")

def insert_database_genomeIDs(database, genomeIDs):
    """
    22.6.24
        Args: 
            Database    Name of the database to be appended
            GenomeIDs   List of genomeIDs to be added to the Genomes table
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute("""PRAGMA synchronous = OFF;""")
        cur.execute("""PRAGMA journal_mode = OFF;""")        
        # Prepare a list of tuples, each containing one genomeID
        genomeID_tuples = [(genomeID,) for genomeID in genomeIDs]
        # Use executemany to insert all genomeIDs in a single batch
        cur.executemany('''INSERT OR IGNORE INTO Genomes (genomeID) VALUES (?)''', genomeID_tuples)
        con.commit()
    con.close()
    return    
        
def insert_database_proteins(database, protein_dict):
    """
        Inserts for concated glob hmmsearches. GenomeId must be defined within the protein object
        Args:
            protein_dict    dictionary with protein objects
    """
def insert_database_proteins(database, protein_dict):
    """
    Inserts for concated glob hmmsearches. GenomeId must be defined within the protein object
    Args:
        protein_dict    dictionary with protein objects
    """
    try:
        with sqlite3.connect(database) as con:
            # Inside 'with con', if an exception occurs, it will rollback automatically
            cur = con.cursor()
            cur.execute("""PRAGMA foreign_keys = ON;""")
            cur.execute("""PRAGMA synchronous = OFF;""")
            cur.execute("""PRAGMA journal_mode = OFF;""")
            
            protein_list = sorted(protein_dict.values(), key=lambda x: (x.gene_contig, x.gene_start))
            protein_records = []
            domain_records = []

            for protein in protein_list:
                genomeID = protein.genomeID
                proteinID = f"{genomeID}-{protein.proteinID}"
                domains = protein.domains
                if not genomeID or not proteinID:
                    print(f"Warning error annotated protein {protein.proteinID}")
                    continue

                # Prepare the protein record
                protein_record = (
                    proteinID, genomeID, protein.gene_locustag, protein.gene_contig,
                    protein.gene_start, protein.gene_end, protein.gene_strand,
                    len(domains), protein.get_sequence()
                )
                protein_records.append(protein_record)

                # Prepare domain records
                for domain in domains.values():
                    domain_record = (
                        proteinID, domain.HMM, domain.start, domain.end, domain.score
                    )
                    domain_records.append(domain_record)

            # Batch insert for proteins
            cur.executemany('''INSERT OR IGNORE INTO Proteins
                (proteinID, genomeID, locustag, contig, start, end, strand, dom_count, sequence)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)''', protein_records)
            
            # Batch insert for domains
            cur.executemany('''INSERT OR IGNORE INTO Domains
                (proteinID, domain, domStart, domEnd, score)
                VALUES (?, ?, ?, ?, ?)''', domain_records)

            # No need to manually call con.commit() — 'with' will handle it

    except Exception as e:
        error_message = f"\nError occurred: {str(e)}"
        traceback_details = traceback.format_exc()
        print(f"\tWARNING: Due to an error - {error_message}")
        print(f"\tTraceback details:\n{traceback_details}")

    return
    
    

def insert_database_clusters(database, cluster_dict):
    """
    22.06.2024
        Args:
            database - path to the database file
            cluster_dict - dictionary with cluster objects
    """
    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            cur.execute("""PRAGMA foreign_keys = ON;""")
            cur.execute("""PRAGMA synchronous = OFF;""")
            cur.execute("""PRAGMA journal_mode = OFF;""") 

            cluster_inserts = []
            protein_updates = []
            keyword_inserts = []

            for cluster_index, cluster in cluster_dict.items():
                clusterID = cluster.get_clusterID()
                if clusterID is None:
                    continue

                cluster_inserts.append((clusterID, cluster.genomeID))
                
                for proteinID in cluster.get_genes():
                    full_proteinID = f"{cluster.genomeID}-{proteinID}"
                    protein_updates.append((clusterID, full_proteinID))

                for Keyword in cluster.get_keywords():
                    keyword_inserts.append((clusterID, Keyword.get_keyword(), Keyword.get_completeness(), Keyword.get_csb()))

            # Batch insert clusters
            cur.executemany("""INSERT OR IGNORE INTO Clusters (clusterID, genomeID) VALUES (?, ?)""", cluster_inserts)

            # Batch update proteins
            cur.executemany("""UPDATE Proteins SET clusterID = ? WHERE proteinID = ?""", protein_updates)

            # Batch insert or replace keywords
            cur.executemany("""INSERT OR REPLACE INTO Keywords (clusterID, keyword, completeness, collinearity) VALUES (?, ?, ?, ?)""", keyword_inserts)
    except Exception as e:
        error_message = f"\nError occurred: {str(e)}"
        traceback_details = traceback.format_exc()
        print(f"\tWARNING: Due to an error - {error_message}")
        print(f"\tTraceback details:\n{traceback_details}")
 

    return

  
def insert_taxonomy_data(database, taxonomy_file):
    """
    Insert taxonomy data from a file into the Genomes table of the database.

    Args:
        database: Path to the SQLite database file.
        taxonomy_file: Path to the parsed taxonomy file (tab-separated).
    """
    try:
        # Check if the database file path is valid
        if not os.path.exists(os.path.dirname(database)):
            raise FileNotFoundError(f"Directory for database does not exist: {os.path.dirname(database)}")

        # Connect to the SQLite database
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            
            # Read the taxonomy file and insert data into the Genomes table
            with open(taxonomy_file, 'r') as file:
                # Skip the header
                next(file)
                
                # Read each line in the file
                for line in file:
                    fields = line.strip().split('\t')
                   
                    if len(fields) < 8:
                        print(f"Skipping line due to insufficient columns: {line}")
                        continue                
                   
                    # Prepare the data for insertion
                    genome_id = fields[0]
                    superkingdom = fields[1]
                    clade = ""
                    phylum = fields[2]
                    class_ = fields[3]
                    order = fields[4]
                    family = fields[5]
                    genus = fields[6]
                    species = fields[7]

                    # Insert the data into the Genomes table
                    cur.execute('''
                        INSERT INTO Genomes (genomeID, Superkingdom, Clade, Phylum, Class, Ordnung, Family, Genus, Species)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ON CONFLICT(genomeID) DO UPDATE SET
                        Superkingdom=excluded.Superkingdom,
                        Clade=excluded.Clade,
                        Phylum=excluded.Phylum,
                        Class=excluded.Class,
                        Ordnung=excluded.Ordnung,
                        Family=excluded.Family,
                        Genus=excluded.Genus,
                        Species=excluded.Species
                    ''', (genome_id, superkingdom, clade, phylum, class_, order, family, genus, species))
            
            # Commit the transaction
            con.commit()
    
    except sqlite3.OperationalError as e:
        print(f"SQLite Operational Error: {e}")
        print(f"Database path: {database}")
        print(f"Ensure the file exists and you have read/write permissions.")
    except sqlite3.Error as e:
        print(f"SQLite Error: {e}")
    except FileNotFoundError as e:
        print(f"File Not Found Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        
        
        
##############################################################
########## Alter information from database routines ##########
##############################################################

def update_domain(database, protein_diction, old_tag, new_tag):
    """
    18.11.22
        Args:
           database     Name of the database to be worked on
           protein_diction Dictionary of protein IDs to be updated
           old_tag     The old domain tag that needs to be updated
           new_tag     The new domain tag that will replace the old tag
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        query = """UPDATE Domains SET domain = ? WHERE proteinID = ? AND domain = ?"""
        
        # Prepare the data for bulk update
        data = [(new_tag, proteinID, old_tag) for proteinID in protein_diction]
        
        # Use executemany for batch processing
        cur.executemany(query, data)
        
        # Commit the transaction (optional, as it is done automatically with 'with' context)
        con.commit()

    return

##############################################################
########## Fetch information from database routines ##########
##############################################################

def update_keywords(database, keyword_dict, batch_size=400):
    """
    Update keywords in the database in batches.
    
    Args:
        database: Name of the database to be worked on.
        keyword_dict: Dictionary containing the new keywords and associated clusterIDs.
        batch_size: Number of (clusterID, keyword) pairs per batch insert.
    """
    if keyword_dict:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            query = """INSERT INTO Keywords (clusterID, keyword) VALUES (?, ?)"""
            
            inserts = []
            for new_keyword, clusterIDs in keyword_dict.items():
                for clusterID in clusterIDs:
                    inserts.append((clusterID, new_keyword))
            
            # In Batches verarbeiten
            for i in range(0, len(inserts), batch_size):
                batch = inserts[i:i+batch_size]
                cur.executemany(query, batch)

    return



def delete_keywords_from_csb(database, options):
    """
    Remove keywords from the database that match the pattern options.csb_name_prefix + a number + options.csb_name_suffix.

    Args:
        database: Name of the database to be worked on.
        options: An object containing csb_name_prefix and csb_name_suffix attributes.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        
        # Construct the pattern
        pattern = f"{options.csb_name_prefix}%{options.csb_name_suffix}"
        
        # SQL query to delete matching keywords
        delete_query = "DELETE FROM Keywords WHERE keyword LIKE ?"
        
        try:
            cur.execute(delete_query, (pattern,))
            con.commit()
        except sqlite3.Error as e:
            print("SQLite error:", e)
            raise

    return



def fetch_genomeIDs(database):
    """
    Fetches the distinct genomeIDs from the Genomes table in the SQLite database.
    
    Args:
        database: Path to the SQLite database.
    
    Returns:
        A set of distinct genomeIDs.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("SELECT DISTINCT genomeID FROM Genomes")
        
        # Use set comprehension to create the set of genomeIDs
        genomeIDs = {row[0] for row in cur.fetchall()}
    
    return genomeIDs


def clean_database_locks(database_path, wait_seconds=10):
    """
    Ensure the database is unlocked and all pending writes are flushed.

    Args:
        database_path: Path to the SQLite .db file
        wait_seconds: How long to wait if database is busy (default 10 seconds)
    """
    wal_path = database_path + "-wal"
    shm_path = database_path + "-shm"

    print(f"Checking database for locks: {database_path}")

    # Step 1: Check if WAL and SHM files exist
    wal_exists = os.path.exists(wal_path)
    shm_exists = os.path.exists(shm_path)

    if wal_exists or shm_exists:
        print(f"Found WAL/SHM files. WAL: {wal_exists}, SHM: {shm_exists}")
    else:
        print("No WAL/SHM files found. Database seems clean.")

    # Step 2: Try to connect and perform a WAL checkpoint
    print("Attempting to checkpoint database")
    success = False
    start_time = time.time()

    while not success and (time.time() - start_time) < wait_seconds:
        try:
            with sqlite3.connect(database_path, timeout=5) as con:
                # Force complete WAL checkpoint
                con.execute("PRAGMA wal_checkpoint(FULL);")
                con.commit()
                success = True
                print("Checkpoint successful. Database flushed.")
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e):
                print("Database is locked. Waiting a bit...")
                time.sleep(1)
            else:
                print(f"Unexpected SQLite error: {e}")
                raise

    if not success:
        raise RuntimeError(f"Could not unlock the database after {wait_seconds} seconds.")

    # Step 3: Check again if WAL and SHM still exist
    wal_exists = os.path.exists(wal_path)
    shm_exists = os.path.exists(shm_path)

    if not wal_exists and not shm_exists:
        print("WAL/SHM files cleaned up successfully.")
    else:
        print("WAL/SHM files still exist. Database might have unclean shutdown earlier.")
        # Optional: Force delete them if you are absolutely sure
        # os.remove(wal_path)
        # os.remove(shm_path)

    print("Database cleaning routine finished.")


