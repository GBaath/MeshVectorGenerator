# MeshVectorGenerator
*A solution for generating points within the bounds of scaled and rotated meshes in Unreal*
The Entire cpp file is located [**here**](https://github.com/GBaath/UnrealFunctionLibraries/tree/main/.cpp) 

*Interpolating points witihn is currently w.i.p. some class members are unused or unimplemented*

Summary: Read all the vertecies, and ranomize points within sphere or box, check rayintersection with the transformed input mesh, optimizing the generation space overtime.
Big & bulky mesh = faster, generating bounding volume not meant for realtime. 

---

# Generate Bounding Volume
**Generate vectors within limits of box or sphere depending on best fit**

**Mesh is quite uniform - boundingsphere is most optimal, yellow = max Radius, red = avarage**
<img src="Radiuses.png" width="50%"/>
<img src="OrientedBox.png" width="50%"/>


<details>
<summary>GetDataFromVertecies</summary>

 ```cpp
//call this before generating point
void UVolumeInterpolator::Init(UStaticMeshComponent* Mesh, const FTransform WT)
{
    FJsonSerializableArrayInt Triangles;
    TArray<FVector> Normals;

    //save vertex data for point in polygon alg
    GetDataFromVerticies(Mesh, WT, Vertecies, Triangles, Normals, AvarageRadius, LargestRadius, BoxvolumeIsSmallest);
}


void UVolumeInterpolator::GetDataFromVerticies(UStaticMeshComponent* Mesh, const FTransform WT, TArray<FVector>& Verts, FJsonSerializableArrayInt& TriangleIndexes, TArray<FVector>& Normals,
    float& OutAvgRadiusToCenter, float& OutMaxRadiusToCenter, bool& bBoxIsSmallestBoundingShape) {

    OutMaxRadiusToCenter = 0;
    OutAvgRadiusToCenter = 0;

    TArray<FVector2D> UV;
    TArray<FProcMeshTangent> Tangents;
    UKismetProceduralMeshLibrary::GetSectionFromStaticMesh(Mesh->GetStaticMesh(), 0, 0, Verts, TriangleIndexes, Normals, UV, Tangents);
    float MagnitudeValues = 0;
    float Radius = 0;

    FOrientedBox Box;


    //REPLACE FUNCTION WITH VARIABLES SET FROM START
    //verts are localspace
    UVolumeInterpolator::GetSmallestBoundingShape(WT, Verts, Radius ,bBoxIsSmallestBoundingShape);
    UE_LOG(LogTemp, Warning, TEXT("%d"), Radius);

    //Randomize point using radius or box depening on smallest error margin 
    if (bBoxIsSmallestBoundingShape) {
        return;
    }
    else {
        //find avarageradius
        for (FVector& vert : Verts) {

            UCommonFunctions::TransformVector(vert, WT);
            //UKismetSystemLibrary::DrawDebugPoint(Mesh, vert, 3.f, FLinearColor::Red, 5.f

            float Magnitude = (vert - WT.GetLocation()).Length();
            MagnitudeValues += Magnitude;
        }
        if (Verts.Num() <= 0) {
            UE_LOG(LogTemp, Warning, TEXT("Verticies Array is 0"));
            return;
        }

        OutAvgRadiusToCenter = MagnitudeValues / Verts.Num();
        OutAvgRadiusToCenter;
        OutMaxRadiusToCenter = Radius;
        OutMaxRadiusToCenter;
    }
    return;
}


```
</details>
<details>
<summary>GetSmallestBoundingShape</summary>

 ```cpp
//this should not be recalculated during runtime
void UVolumeInterpolator::GetSmallestBoundingShape(const FTransform WT, TArray<FVector> Verts, float& OutRadius , bool& BoxIsSmallest) {


    UE::Geometry::TMinVolumeBox3<float> box;

    Vertecies = Verts;


    TFunctionRef<FVector3f(int32)> GetVertexInWorldSpace =
        [Verts,WT](int32 i) {

        FVector V  = Verts[FMath::Clamp(i, 0, Verts.Num() - 1)];
        UCommonFunctions::TransformVector(V, WT);

        return (FVector3f)V;
        };

    //get smallest boundingbox
    UE::Geometry::FOrientedBox3f OBox;

    //Solve is somewhat unreliable
    int i = 0;
    int maxTries = 10;
    do {
        box.Solve(Verts.Num(), GetVertexInWorldSpace, true);
        i++;
    } while (i < maxTries || !box.IsSolutionAvailable());

    if (box.IsSolutionAvailable()) {
        box.GetResult(OBox);
    }
    else
    {
        UE_LOG(LogTemp, Warning, TEXT("Failed To Solve BoxVolume"));

        BoxIsSmallest = false;
        //return SphereData
        //largest radius is smallest boundingbox
        float LRadius = 0.f;
        for (FVector& vert : Verts) {
            FVector V = vert;
            UCommonFunctions::TransformVector(V, WT);
            //UCommonFunctions::TransformVector(V, WorldTransform);
            UKismetSystemLibrary::DrawDebugPoint(GetWorld(), V, 5.f, FLinearColor::White, 1.f);
            LRadius = FMath::Max((V-WT.GetLocation()).Length(), LRadius);
        }

       BoxvolumeIsSmallest = BoxIsSmallest;
       OutRadius = LRadius;
       LargestRadius = LRadius;


        return;
    }
    //largest radius is smallest bounding
    float LRadius = 0.f;
    for (FVector& vert : Verts) {
        FVector V = vert;
        UCommonFunctions::TransformVector(V, WT);
        //UCommonFunctions::TransformVector(V, WorldTransform);
        UKismetSystemLibrary::DrawDebugPoint(GetWorld(), V, 5.f, FLinearColor::White, 1.f);
        LRadius = FMath::Max((V-WT.GetLocation()).Length(), LRadius);
    }

    LargestRadius = LRadius;

    float boxVol = OBox.Volume();
    float radVol = ((4 / 3) * UKismetMathLibrary::GetPI() * FMath::Cube(LRadius));
    //compare volumes
    if (OBox.Volume() > ((4 / 3) * UKismetMathLibrary::GetPI() * FMath::Cube(LRadius))) {
        BoxIsSmallest = false;
        OutRadius = LRadius;

        BoxvolumeIsSmallest = BoxIsSmallest;

        return;
    }

    BoxIsSmallest = true;


    OrientedBox = ConvertFrom3f(OBox);
    BoxvolumeIsSmallest = BoxIsSmallest;

    return;
}

```
</details>

# Cross reference against mesh


# Optimise size

