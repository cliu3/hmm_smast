function db = flipdb(db,dir)

if isfield(db,'tidalro')
    db.tidro = db.tidalro;
    clear db.tidalro
end
if isfield(db,'rough')
    db.bathro = db.rough;
    clear db.rough
end

switch dir
    case 'lat'
        db.lat    = db.lat(end:-1:1,:);
        db.long   = db.long(end:-1:1,:);
        db.depth  = db.depth(end:-1:1,:);
        db.land   = db.land(end:-1:1,:);
        db.bathro = db.bathro(end:-1:1,:);
        db.tidro  = db.tidro(end:-1:1,:);
        db.amp    = db.amp(end:-1:1,:,:);
        db.phase  = db.phase(end:-1:1,:,:);
    case 'long'
        db.lat    = db.lat(:,end:-1:1);
        db.long   = db.long(:,end:-1:1);
        db.depth  = db.depth(:,end:-1:1);
        db.land   = db.land(:,end:-1:1);
        db.bathro = db.bathro(:,end:-1:1);
        db.tidro  = db.tidro(:,end:-1:1);
        db.amp    = db.amp(:,end:-1:1,:);
        db.phase  = db.phase(:,end:-1:1,:);
end