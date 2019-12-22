import boto
import boto3
import pandas as pd

#twinz = the bucket for the msu twins
bucket_name = 'twinz01'
s3 = boto3.resource('s3')
bucket = s3.Bucket(bucket_name)
images = [img.key for img in bucket.objects.all()]
client = boto3.client('rekognition')

results_long = []

#'img1' or 'img2' refers to the img in the 'images' list, which is the collection images from the bucket;
#so, this script will generate a comparison between every image in the bucket against every other image.
#For example, if you have 10 images in a bucket, you will generate 10 x 10 = 100 comparisons,
#including the image compared to itself (negative control).

#for img1 in images:
#    for img2 in images:
for x in range(0,len(images)-1):
    for y in range( (x+1),len(images)):
    
        image1=images[x]
        image2=images[y]
        
    
        response = client.compare_faces(
            SourceImage={
                'S3Object': {
                    'Bucket': bucket_name,
                    'Name': image1
                }
            },
            TargetImage={
                'S3Object': {
                    'Bucket': bucket_name,
                    'Name': image2
                }
        },
        SimilarityThreshold=0 
        )
        if 'FaceMatches' in response:
            for f, face in enumerate(response['FaceMatches']):    
                results_long.append(
                {'target': image1, 
                'source': image2,
                'Similarity': face['Similarity']})
    
            
img_df_long = pd.DataFrame(results_long, columns=['target','source','Similarity'])
csv = "schoeller_cf_nonredundant.csv"
img_df_long.to_csv(csv)
