class Model 
{
public:
    Model ();
    virtual ~Model();
	
private:
    float* mVertices;
	float* mVerticeNormals;
	unsigned int* mFaces;
	
	unsigned int mNumVertices;
	unsigned int mNumFaces;	
};
