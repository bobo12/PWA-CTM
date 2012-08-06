function data = collectDataFromDB(dataPSName, dataPS, columnTypes, ...
                                  setData, setTypes)

    readDB = core.DatabaseReader('localhost', 'pems');
    readDB.psCreate(dataPSName, dataPS);
    numColumns = length(columnTypes);
    numSetColumns = length(setTypes);

    
    for iColumn = 1:numSetColumns
        setter = setterColumnType(readDB, setTypes{iColumn});
        setter(dataPSName, iColumn, setData{iColumn});
    end
    
    readDB.psQuery(dataPSName);
    
    numRows = readDB.getFetchSize;
    columnNames = readDB.psRSColumnNames(dataPSName);
    
    data = cell(numRows, numColumns);

    iRow = 1;
    while (readDB.psRSNext(dataPSName))
        for iColumn = 1:numColumns
            parser = getterColumnType(readDB, columnTypes{iColumn});
            data{iRow, iColumn} = parser(dataPSName, ...
                                         columnNames(iColumn));
        end
        iRow = iRow + 1;
    end
    readDB.psDestroy(dataPSName);
    readDB.close();
end

function dbMethod = setterColumnType(dbReader, setTypes)
    switch(setTypes)
      case 'int'
        dbMethod = @dbReader.psSetInteger;
      case 'boolean'
        dbMethod = @dbReader.psSetBoolean;
      case 'real'
        dbMethod = @dbReader.psSetReal;
      case 'date'
        dbMethod = @dbReader.psSetTimestamp;
      case 'long'
        dbMethod = @dbReader.psSetBigIntl
      case 'short'
        dbMethod = @dbReader.psSetSmallInt;
      case 'aint'
        dbMethod = @dbReader.psSetArrayInteger;
      case 'areal'
        dbMethod = @dbReader.psSetArrayReal;
      case 'ashort'
        dbMethod = @dbReader.psSetArraySmallInt;
      case 'string'
        dbMethod = @dbReader.psSetVarChar;
      case 'geom'
        dbMethod = @dbReader.psSetPostGISGeometry;
      case 'multiline'
        dbMethod = @dbReader.psSetPostGISLineString;
    end

end

function dbMethod = getterColumnType(dbReader, columnType)
    switch(columnType)
      case 'int'
        dbMethod = @dbReader.psRSGetInteger;
      case 'boolean'
        dbMethod = @dbReader.psRSGetBoolean;
      case 'real'
        dbMethod = @dbReader.psRSGetReal;
      case 'date'
        dbMethod = @dbReader.psRSGetTimestamp;
      case 'long'
        dbMethod = @dbReader.psRSGetBigInt;
      case 'short'
        dbMethod = @dbReader.psRSGetSmallInt;
      case 'aint'
        dbMethod = @dbReader.psRSGetArrayInteger;
      case 'areal'
        dbMethod = @dbReader.psRSGetArrayReal;
      case 'ashort'
        dbMethod = @dbReader.psRSGetArraySmallInt;
      case 'string'
        dbMethod = @dbReader.psRSGetVarChar;
      case 'geom'
        dbMethod = @dbReader.psRSGetPostGISGeometry;
      case 'multiline'
        dbMethod = @dbReader.psRSGetPostGISLineString;
    end
end