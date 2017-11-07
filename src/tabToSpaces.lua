local i = 1
local inFile = arg[i]
while arg[i] ~= nil do
    io.input(arg[i])
    local text = io.read("*all")
    local outText = ""
    
    for j = 0, #text do
        if text:sub(j,j) == "\t" then
            outText = outText .. "    "
        else
            outText = outText .. text:sub(j,j)
        end
    end
    io.output(arg[i] .. ".no-tabs")
    io.write(outText)
    i = i + 1
end
