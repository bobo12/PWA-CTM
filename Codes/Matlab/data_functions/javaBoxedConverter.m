function doubleNum = javaBoxedConverter(javaNumber)
    if isjava(javaNumber) && strcmp(...
            javaNumber.getClass.getSuperclass, ...
            'class java.lang.Number')
        doubleNum = javaNumber.doubleValue;
    elseif ~isjava(javaNumber)
        doubleNum = javaNumber;
    else
        error('Not a number!\n');
    end
end