#!/usr/bin/env python3

import sqlite3
import argparse


def list_tables(cur):
    cur.execute("""
        SELECT name
        FROM sqlite_master
        WHERE type = 'table'
        ORDER BY name;
    """)
    return [row[0] for row in cur.fetchall()]


def get_columns(cur, table):
    cur.execute(f'PRAGMA table_info("{table}")')
    return [row[1] for row in cur.fetchall()]


def print_table(cur, table, limit=None):
    columns = get_columns(cur, table)

    print("\n" + "=" * 80)
    print(f"TABLE: {table}")
    print("=" * 80)

    if not columns:
        print("[EMPTY / NO COLUMNS]")
        return

    print("\t".join(columns))

    query = f'SELECT * FROM "{table}"'
    if limit is not None:
        query += f" LIMIT {int(limit)}"

    cur.execute(query)

    rows = cur.fetchall()

    if not rows:
        print("[NO ROWS]")
        return

    for row in rows:
        row = ["NULL" if value is None else str(value) for value in row]
        print("\t".join(row))


def dump_database(database, limit=None):
    con = sqlite3.connect(database)
    cur = con.cursor()

    tables = list_tables(cur)

    print(f"[INFO] Database: {database}")
    print(f"[INFO] Number of tables: {len(tables)}")

    for table in tables:
        print_table(cur, table, limit=limit)

    con.close()


def main():
    parser = argparse.ArgumentParser(
        description="Print all tables and their content from a SQLite database."
    )

    parser.add_argument(
        "-db",
        "--database",
        required=True,
        help="Path to SQLite database."
    )

    parser.add_argument(
        "-l",
        "--limit",
        type=int,
        default=None,
        help="Optional maximum number of rows printed per table."
    )

    args = parser.parse_args()

    dump_database(args.database, limit=args.limit)


if __name__ == "__main__":
    main()