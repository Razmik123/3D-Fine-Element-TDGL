function model = interfaceModel(regionIds, orderParameterType, options)
%INTERFACEMODEL Declare the order-parameter law for a material interface.

arguments
    regionIds (1,2) double {mustBeInteger,mustBePositive}
    orderParameterType (1,1) string {mustBeMember(orderParameterType, ...
        ["continuous","de-gennes","no-order-parameter"])}
    options.GammaB (1,1) double {mustBeNonnegative} = 0
end

if regionIds(1) == regionIds(2)
    error('tdgl:problem:InvalidInterfaceRegions', ...
        'An interface must connect two different material regions.');
end
model = struct( ...
    'regionIds',sort(regionIds), ...
    'orderParameterType',orderParameterType, ...
    'gammaB',options.GammaB);
end
