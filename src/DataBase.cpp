#include "SOM.hpp"

DataBase::DataBase()
{
	;
}

DataBase::~DataBase()
{
	;
}

int DataBase::open(char *fileName)
{
	int rc = sqlite3_open_v2(fileName, &db, SQLITE_OPEN_READONLY, NULL);
	
	if( verbose )
		std::cout << "Opening database: " << fileName << "\n";
	
	// If error
    if (rc != SQLITE_OK) {
        
        std::cout << "Cannot open database: " << sqlite3_errmsg(db) << "\n";
        sqlite3_close(db);
        
        return (FALSE);
    }
	return TRUE;
}

void DataBase::close()
{
	
}

int DataBase::getRecord(char *buff[], int row, std::string icanName)
{
	return 0;
}

int DataBase::startTransaction()
{
	char *err_msg;
	int retValue{};
	if( retValue = sqlite3_exec(db, "BEGIN;", NULL, NULL, &err_msg) )
	{
		fprintf(stderr, "Failed to BEGIN transaction: %s\n", sqlite3_errmsg(db));
	}
	return retValue;
}

int DataBase::endTransaction()
{
	char *err_msg;
	int retValue{};
	if( retValue = sqlite3_exec(db, "END;", NULL, NULL, &err_msg) )
	{
		fprintf(stderr, "Failed to END transaction: %s\n", sqlite3_errmsg(db));
	}
	return retValue;
}

int DataBase::getElement(char *buff, int row, std::string icanName)
{
	int rc;
	char sql[200];
	char *err_msg;
	
	sprintf(sql, "SELECT %s FROM Ican WHERE Id = '%d';", icanName.c_str(), row);
	
	// Exekvera SQL statement
	rc = sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if (rc != SQLITE_OK ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
		return(-1);
    }
	
	return(0);
}

int saveElement(void *buff, int argc, char **argv, char **azColName)
{
	char *ptr = (char*)buff;
	
	// If value string is not NULL
	if( argv[0] )
	{
		strcpy(ptr, argv[0]);
		//std::cout << azColName[0] << ":\t" << ptr << "\n";
	}
	else
	{
		strcpy(ptr, "NULL");
	}
	return 0;
}

int DataBase::rows()
{
	int rc;
	char *err_msg;
	char sql[] = "SELECT COUNT(*) FROM Ican;";
	char buff[50];
	unsigned long value;
	
	// Exekvera SQL statement
	rc = sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if (rc != SQLITE_OK ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
		return(-1);
    }
	
	value = strtoul(buff, NULL, 10);
	
	return((int)value);
}

int DataBase::minId()
{
	int rc;
	char *err_msg;
	char sql[] = "SELECT MIN(Id) FROM Ican;";
	char buff[50];
	unsigned long value;
	
	// Exekvera SQL statement
	rc = sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if (rc != SQLITE_OK ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
		return(-1);
    }
	
	value = strtoul(buff, NULL, 10);
	
	return((int)value);
}

int DataBase::maxId()
{
	int rc;
	char *err_msg;
	char sql[] = "SELECT MAX(Id) FROM Ican;";
	char buff[50];
	unsigned long value;
	
	// Exekvera SQL statement
	rc = sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if (rc != SQLITE_OK ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
		return(-1);
    }
	
	value = strtoul(buff, NULL, 10);
	
	return((int)value);
}

int DataBase::doesExist(int row)
{
	int rc;
	char *err_msg;
	char sql[200];
	char buff[50];
	unsigned long value;
	
	sprintf(sql, "SELECT COUNT(*) FROM Ican WHERE Id = \'%d\';", row);
	
	// Exekvera SQL statement
	rc = sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if (rc != SQLITE_OK ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
		return(-1);
    }
	
	value = strtoul(buff, NULL, 10);
	
	//if( !value ) printf("\nDoes not exist\n");
	
	return((int)value);
}

int DataBase::getMax(char *buff, std::string icanName)
{
	int rc;
	char sql[200];
	char *err_msg;
	
	sprintf(sql, "SELECT MAX(CAST(%s AS DECIMAL)) FROM Ican;", icanName.c_str());
	
	// Exekvera SQL statement
	rc = sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if (rc != SQLITE_OK ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
		return(-1);
    }
	
	return(0);
}