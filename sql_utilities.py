try:
    import mysql.connector as db_connector
except:
    print "MySQLdb import error, make sure that it is properly installed!"


def example_chembl_query():
    dbhost = None 
    dbuser = "emre" 
    dbname = "chembl_23"
    cursor, db = sql_utilities.get_db_cursor(dbhost, dbname, dbuser, dbpass=None)
    tables = ["compound_records AS C", "molecule_dictionary AS M"]
    columns = ["M.chembl_id"] 
    fixed_conditions = [["C.compound_name", "=", "donepezil"]]
    join_conditions = [["C.molregno", "=", "M.molregno"]]
    sql = sql_utilities.get_select_sql_query(tables, columns = columns, fixed_conditions = fixed_conditions, join_conditions = join_conditions, group_conditions=None, distinct_columns=True)
    print sql
    results = sql_utilities.get_values_from_database(cursor, sql)
    if results:
	for row in results:
	    print row
    sql_utilities.close_db_cursor(cursor, db)
    return


def get_db_cursor(dbhost, dbname, dbuser, dbpass=None):
    db = db_connector.connect(host = dbhost, user = dbuser, passwd = dbpass) # port = dbport, unix_socket= dbsocket)
    cursor = db.cursor()
    cursor.execute("USE %s" % dbname)
    #cursor.execute("""Select LAST_INSERT_ID()""")
    return cursor, db


def close_db_cursor(cursor, db):
    cursor.close()
    db.close()
    return


def get_values_from_database(cursor, sql): 
    cursor.execute(sql)
    answer = cursor.fetchall()
    return answer


def get_select_sql_query(tables, columns=None, fixed_conditions=None, join_conditions=None, group_conditions=None, distinct_columns=False):
    """
    Generates a general select sql statement
    "tables" is a list or tuple of tables where the value/s must be searched. If the elements of the list or tuple are tuples of length 2, the format taken will be the following:
		  (table_name or table_object, alias to the table)
    "columns" is a list or tuple of the columns searched (columns must be preceeded by the table where they are searched). If it is None, all values will be selected
    "fixed_conditions" is a list or tuple of tuples with the following format: (column,type,restriction_value)
    "join_conditions" is a list or tuple of tuples with the following format: (column,type,column) to restrict the selection to the joint
    "type" can be "=",">","<",...
    "group_conditions" is a list of columns where it must be grouped 
    It returns the sql query.
    Method adapted from BIANA
    """
    if( fixed_conditions is None or not len(fixed_conditions) ):
	fixed_conditions_sql = ""
    else:
	fixed_conditions_sql = " AND ".join(["%s %s \"%s\"" %(x[0],x[1],x[2]) for x in fixed_conditions])

    if( join_conditions is None or not len(join_conditions) ):
	join_conditions_sql = ""
    else:
	join_conditions_sql = " AND ".join(["%s %s %s" %(x[0],x[1],x[2]) for x in join_conditions])
	if fixed_conditions_sql != "":
	    join_conditions_sql = " AND %s" %(join_conditions_sql)

    if( join_conditions or fixed_conditions ):
	where_sql = " WHERE "
    else:
	where_sql = ""

    if columns is None:
	columns_sql = "*"
    else:
	columns_list = []
	for current_column in columns:
	    if( isinstance(current_column, tuple) ):
		columns_list.append("%s AS %s" %(current_column))
	    else:
		columns_list.append(current_column)
	columns_sql = ",".join(columns_list)

    # tranform table objects to table name strings
    tables_list = []
    for actual_table in tables:
	if( isinstance(actual_table,tuple) ):
	    tables_list.append("%s AS %s " %(actual_table[0],actual_table[1]) )
	else:
	    tables_list.append("%s" %actual_table)
			    
    if group_conditions is not None and len(group_conditions)>0:
	group_conditions_sql = "GROUP BY %s" %(",".join(group_conditions))
    else:
	group_conditions_sql = ""

    if distinct_columns:
	distinct_str = "DISTINCT"
    else:
	distinct_str = ""
	
    return """SELECT %s %s FROM %s %s %s %s %s""" %(distinct_str,
						    columns_sql,
						    ",".join(tables_list),
						    where_sql,
						    fixed_conditions_sql,
						    join_conditions_sql,
						    group_conditions_sql)
 
