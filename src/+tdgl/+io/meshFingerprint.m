function fingerprint = meshFingerprint(mesh)
%MESHFINGERPRINT Stable SHA-256 identity for mesh geometry and regions.

digest = java.security.MessageDigest.getInstance('SHA-256');
updateDigest(digest,double(mesh.nodes(:)));
updateDigest(digest,uint32(mesh.cells(:)));
updateDigest(digest,int32(mesh.regionIds(:)));
raw = typecast(digest.digest(),'uint8');
fingerprint = lower(string(reshape(dec2hex(raw,2).',1,[])));
end

function updateDigest(digest,value)
bytes = typecast(value(:),'uint8');
digest.update(typecast(bytes,'int8'));
end
