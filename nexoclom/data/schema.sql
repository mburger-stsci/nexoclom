CREATE TABLE geometry_with_time (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    Planet SSObject NOT NULL,
    StartPoint SSObject NOT NULL,
    objects SSObject[],
    starttime TIMESTAMP NOT NULL
);

CREATE TABLE geometry_without_time (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    Planet SSObject NOT NULL,
    StartPoint SSObject NOT NULL,
    objects SSObject ARRAY NOT NULL,
    phi DOUBLE PRECISION ARRAY,
    subsolarpt DOUBLE PRECISION[2] NOT NULL,
    TAA DOUBLE PRECISION NOT NULL,
    check (cardinality(objects)-1 = cardinality(phi)),
    check (TAA >= 0 and TAA <= 6.284)
    -- can add constraint on subsolarpoint and phi
);

CREATE TABLE surface_int_constant (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    stickcoef DOUBLE PRECISION NOT NULL,
    accomfactor DOUBLE PRECISION,
    check (stickcoef >= 0 and stickcoef <= 1),
    check (accomfactor >= 0 and accomfactor <= 1)
);

CREATE TABLE surface_int_map (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    mapfile TEXT NOT NULL,
    accomfactor DOUBLE PRECISION NOT NULL,
    check (accomfactor >= 0 and accomfactor <= 1)
);

CREATE TABLE surface_int_tempdependent (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    accomfactor DOUBLE PRECISION NOT NULL,
    a DOUBLE PRECISION[3] NOT NULL,
    check (accomfactor >= 0 and accomfactor <= 1)
);

CREATE TABLE forces (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    gravity BOOLEAN NOT NULL,
    radpres BOOLEAN NOT NULL
);

CREATE TABLE spatdist_uniform (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    exobase DOUBLE PRECISION NOT NULL,
    longitude DOUBLE PRECISION[2] NOT NULL,
    latitude DOUBLE PRECISION[2] NOT NULL,
    check(exobase >= 1)
);

CREATE TABLE spatdist_surfmap (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    exobase DOUBLE PRECISION NOT NULL,
    mapfile TEXT NOT NULL,
    coordsystem TEXT NOT NULL,
    check(exobase >= 1)
);

CREATE TABLE spatdist_spot (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    exobase DOUBLE PRECISION NOT NULL,
    longitude DOUBLE PRECISION NOT NULL,
    latitude DOUBLE PRECISION NOT NULL,
    sigma DOUBLE PRECISION NOT NULL,
    check(exobase >= 1)
);

CREATE TABLE speeddist_gaussian (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    vprob DOUBLE PRECISION NOT NULL,
    sigma DOUBLE PRECISION NOT NULL,
    check (sigma >= 0.)
);

CREATE TABLE speeddist_maxwellian (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    temperature DOUBLE PRECISION NOT NULL,
    check (temperature >= 0)
);

CREATE TABLE speeddist_sputtering (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    alpha DOUBLE PRECISION NOT NULL,
    beta DOUBLE PRECISION NOT NULL,
    U DOUBLE PRECISION NOT NULL
);

CREATE TABLE speeddist_flat (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    vprob DOUBLE PRECISION NOT NULL,
    delv DOUBLE PRECISION NOT NULL,
    check (delv >= 0.)
);

CREATE TABLE angdist_isotropic (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    altitude DOUBLE PRECISION[2] NOT NULL,
    azimuth DOUBLE PRECISION[2] NOT NULL
);

CREATE TABLE options (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    endtime DOUBLE PRECISION NOT NULL,
    species TEXT NOT NULL,
    lifetime DOUBLE PRECISION NOT NULL,
    outer_edge DOUBLE PRECISION NOT NULL,
    step_size DOUBLE PRECISION NOT NULL,
    resolution DOUBLE PRECISION,
    check (endtime > 0),
    check (outer_edge > 0),
    check (step_size >= 0)
);

CREATE TABLE outputfile (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    filename TEXT UNIQUE NOT NULL,
    npackets INT NOT NULL,
    totalsource DOUBLE PRECISION NOT NULL,
    creation_time TIMESTAMP NOT NULL,
    geo_type TEXT NOT NULL,  --"with time" or "without time"
    geo_id INT NOT NULL,
    sint_type TEXT NOT NULL, --"constant", "temperature dependent", "from map"
    sint_id INT NOT NULL,
    force_id INT NOT NULL,
    spatdist_type TEXT NOT NULL,
    spatdist_id INT NOT NULL,
    spddist_type TEXT NOT NULL,
    spddist_id INT NOT NULL,
    angdist_type TEXT NOT NULL,
    angdist_id INT NOT NULL,
    opt_id INT NOT NULL,
    check (npackets > 0),
    check (totalsource > 0)
);

CREATE TABLE modelimages (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    out_idnum INT NOT NULL,
    quantity TEXT NOT NULL,
    origin TEXT NOT NULL,
    dims INT[2] NOT NULL,
    center DOUBLE PRECISION[2] NOT NULL,
    width DOUBLE PRECISION[2] NOT NULL,
    subobslongitude DOUBLE PRECISION NOT NULL,
    subobslatitude DOUBLE PRECISION NOT NULL,
    mechanism TEXT,
    wavelength TEXT,
    filename TEXT NOT NULL)

CREATE TABLE uvvsmodels (
    idnum INT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    out_idnum INT NOT NULL,
    quantity TEXT NOT NULL,
    orbit INT NOT NULL,
    dphi DOUBLE PRECISION NOT NULL,
    mechanism TEXT,
    wavelength TEXT,
    filename TEXT NOT NULL)

DONE
