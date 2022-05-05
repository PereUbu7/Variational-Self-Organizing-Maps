#include "SqliteDataLoader.hpp"

#include <iostream>
#include <cstring>

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

SqliteDataLoader::~SqliteDataLoader()
{
	sqlite3_close(db);
}

bool SqliteDataLoader::open(const char *fileName)
{
	int rc = sqlite3_open_v2(fileName, &db, SQLITE_OPEN_READONLY, NULL);
	
	if( _verbose )
		std::cout << "Opening database: " << fileName << "\n";
	
	// If error
    if (rc != SQLITE_OK) {
        
        std::cout << "Cannot open database: " << sqlite3_errmsg(db) << "\n";
        sqlite3_close(db);
        
        return false;
    }
	return true;
}

void SqliteDataLoader::close()
{
	sqlite3_close(db);
}

int SqliteDataLoader::getRecord(char *buff[], int row, std::string icanName) const
{
	return 0;
}

int SqliteDataLoader::startTransaction()
{
	char *err_msg;
	int retValue{};
	retValue = sqlite3_exec(db, "BEGIN;", NULL, NULL, &err_msg);
	if ( err_msg ) 
	{
		fprintf(stderr, "Failed to BEGIN transaction: Error message:%s\n", err_msg);
		sqlite3_free(err_msg);
		return(-1);
    }
	return retValue;
}

int SqliteDataLoader::endTransaction()
{
	char *err_msg;
	int retValue{};
	retValue = sqlite3_exec(db, "END;", NULL, NULL, &err_msg);
	if ( err_msg ) 
	{
		fprintf(stderr, "Failed to END transaction: Error message:%s\n", err_msg);
		sqlite3_free(err_msg);
		return(-1);
    }
	return retValue;
}

int SqliteDataLoader::getElement(char *buff, int row, std::string icanName)
{
	char sql[200];
	char *err_msg;
	
	sprintf(sql, "SELECT %s FROM Ican WHERE Id = '%d';", icanName.c_str(), row);
	
	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if ( err_msg ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return(-1);
    }
	
	return(0);
}

int SqliteDataLoader::rows()
{
	char *err_msg;
	char sql[] = "SELECT COUNT(*) FROM Ican;";
	char buff[50];
	unsigned long value;
	
	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if ( err_msg ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return(-1);
    }
	
	value = strtoul(buff, NULL, 10);
	
	return((int)value);
}

int SqliteDataLoader::minId()
{
	char *err_msg;
	char sql[] = "SELECT MIN(Id) FROM Ican;";
	char buff[50];
	unsigned long value;
	
	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if ( err_msg ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return(-1);
    }
	
	value = strtoul(buff, NULL, 10);
	
	return((int)value);
}

int SqliteDataLoader::maxId()
{
	char *err_msg;
	char sql[] = "SELECT MAX(Id) FROM Ican;";
	char buff[50];
	unsigned long value;
	
	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if ( err_msg ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return(-1);
    }
	
	value = strtoul(buff, NULL, 10);
	
	return((int)value);
}

int SqliteDataLoader::doesExist(int row)
{
	char *err_msg;
	char sql[200];
	char buff[50];
	unsigned long value;
	
	sprintf(sql, "SELECT COUNT(*) FROM Ican WHERE Id = \'%d\';", row);
	
	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if ( err_msg ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return(-1);
    }
	
	value = strtoul(buff, NULL, 10);
	
	//if( !value ) printf("\nDoes not exist\n");
	
	return((int)value);
}

int SqliteDataLoader::getMax(char *buff, std::string icanName)
{
	char sql[200];
	char *err_msg;
	
	sprintf(sql, "SELECT MAX(CAST(%s AS DECIMAL)) FROM Ican;", icanName.c_str());
	
	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);
    
    if ( err_msg ) 
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return(-1);
    }
	
	return(0);
}