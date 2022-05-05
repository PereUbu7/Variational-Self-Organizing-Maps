#pragma once

#include <string>

class IDataLoader
{
	public:
		virtual ~IDataLoader() = default;
		IDataLoader() = default;
		IDataLoader(const IDataLoader&) = default;
		IDataLoader& operator=(const IDataLoader&) = default;
		IDataLoader(IDataLoader&&) = default;
		IDataLoader& operator=(IDataLoader&&) = default;
		virtual int getElement(char*, int, std::string) = 0;
		virtual int rows() = 0;
		virtual int startTransaction() = 0;
		virtual int endTransaction() = 0;
		virtual int minId() = 0;
		virtual int maxId() = 0;
		virtual int doesExist(int) = 0;
		
};