    Info<< "Reading lowDimensionProperties\n" << endl;

    IOdictionary lowDimensionProperties
    (
        IOobject
        (
            "lowDimensionProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

   Info<< lowDimensionProperties << endl;

   word fieldName=lowDimensionProperties.lookup("lowDimensionField"); 

 Info<< "Reading field selected: " << fieldName<<" field"<<endl<<endl;

    volScalarField C
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
 
 











