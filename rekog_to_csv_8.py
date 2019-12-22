#libraries
import boto
import boto3
import pandas as pd

#the bucket name refers to the name of the database folder, which contains the images
#make sure the name is correct or else you may have the wrong data or generate an error
bucket_name = 'schoeller-ctrl'

#s3 = the database with the buckets
s3 = boto3.resource('s3')

#call upon the bucket with the images that you want to analyze
bucket = s3.Bucket(bucket_name)

#create  variable with all of the images in a list
images = [img.key for img in bucket.objects.all()]

#call upon the facial recognition api
client = boto3.client('rekognition')

#before running the for loop, a list needs to be generated to append the coordinate data
#per image
results_long = []

#the for loop which goes image by image in the 'images' list created on line 16
for img in images:
    response = client.detect_faces(Image={
            'S3Object': {
                'Bucket': bucket_name,
                'Name': img,
            }
        },
        Attributes=['ALL']
    )
    #here, each of the facial landmark coordinates are extracted and normalized to the nose
    #furthermore, they are appended to the 'results_long' list
    if 'FaceDetails' in response:
        for f, face in enumerate(response['FaceDetails']):    
            results_long.append(
            {"ID":img,
            "MouthLeft":[face['Landmarks'][2]['X']-face['Landmarks'][4]['X'],face['Landmarks'][2]['Y']-face['Landmarks'][4]['Y']],
            "MouthRight":[face['Landmarks'][3]['X']-face['Landmarks'][4]['X'],face['Landmarks'][3]['Y']-face['Landmarks'][4]['Y']],
            "LeftEyebrowLeft":[face['Landmarks'][5]['X']-face['Landmarks'][4]['X'],face['Landmarks'][5]['Y']-face['Landmarks'][4]['Y']],
            "LeftEyebrowRight":[face['Landmarks'][6]['X']-face['Landmarks'][4]['X'],face['Landmarks'][6]['Y']-face['Landmarks'][4]['Y']],
            "LeftEyebrowUp":[face['Landmarks'][7]['X']-face['Landmarks'][4]['X'],face['Landmarks'][7]['Y']-face['Landmarks'][4]['Y']],
            "RightEyebrowLeft":[face['Landmarks'][8]['X']-face['Landmarks'][4]['X'],face['Landmarks'][8]['Y']-face['Landmarks'][4]['Y']],
            "RightEyebrowRight":[face['Landmarks'][9]['X']-face['Landmarks'][4]['X'],face['Landmarks'][9]['Y']-face['Landmarks'][4]['Y']],
            "RightEyebrowUp":[face['Landmarks'][10]['X']-face['Landmarks'][4]['X'],face['Landmarks'][10]['Y']-face['Landmarks'][4]['Y']],
            "LeftEyeLeft":[face['Landmarks'][11]['X']-face['Landmarks'][4]['X'],face['Landmarks'][11]['Y']-face['Landmarks'][4]['Y']],
            "LeftEyeRight":[face['Landmarks'][12]['X']-face['Landmarks'][4]['X'],face['Landmarks'][12]['Y']-face['Landmarks'][4]['Y']],
            "LeftEyeUp":[face['Landmarks'][13]['X']-face['Landmarks'][4]['X'],face['Landmarks'][13]['Y']-face['Landmarks'][4]['Y']],
            "LeftEyeDown":[face['Landmarks'][14]['X']-face['Landmarks'][4]['X'],face['Landmarks'][14]['Y']-face['Landmarks'][4]['Y']],
            "RightEyeLeft":[face['Landmarks'][15]['X']-face['Landmarks'][4]['X'],face['Landmarks'][15]['Y']-face['Landmarks'][4]['Y']],
            "RightEyeRight":[face['Landmarks'][16]['X']-face['Landmarks'][4]['X'],face['Landmarks'][16]['Y']-face['Landmarks'][4]['Y']],
            "RightEyeUp":[face['Landmarks'][17]['X']-face['Landmarks'][4]['X'],face['Landmarks'][17]['Y']-face['Landmarks'][4]['Y']],
            "RightEyeDown":[face['Landmarks'][18]['X']-face['Landmarks'][4]['X'],face['Landmarks'][18]['Y']-face['Landmarks'][4]['Y']],
            "NoseLeft":[face['Landmarks'][19]['X']-face['Landmarks'][4]['X'],face['Landmarks'][19]['Y']-face['Landmarks'][4]['Y']],
            "NoseRight":[face['Landmarks'][20]['X']-face['Landmarks'][4]['X'],face['Landmarks'][20]['Y']-face['Landmarks'][4]['Y']],
            "MouthUp":[face['Landmarks'][21]['X']-face['Landmarks'][4]['X'],face['Landmarks'][21]['Y']-face['Landmarks'][4]['Y']],
            "MouthDown":[face['Landmarks'][22]['X']-face['Landmarks'][4]['X'],face['Landmarks'][22]['Y']-face['Landmarks'][4]['Y']],
            "LeftPupil":[face['Landmarks'][23]['X']-face['Landmarks'][4]['X'],face['Landmarks'][23]['Y']-face['Landmarks'][4]['Y']],
            "RightPupil":[face['Landmarks'][24]['X']-face['Landmarks'][4]['X'],face['Landmarks'][24]['Y']-face['Landmarks'][4]['Y']],
            "UpperJawlineLeft":[face['Landmarks'][25]['X']-face['Landmarks'][4]['X'],face['Landmarks'][25]['Y']-face['Landmarks'][4]['Y']],
            "MidJawlineLeft":[face['Landmarks'][26]['X']-face['Landmarks'][4]['X'],face['Landmarks'][26]['Y']-face['Landmarks'][4]['Y']],
            "ChinBottom":[face['Landmarks'][27]['X']-face['Landmarks'][4]['X'],face['Landmarks'][27]['Y']-face['Landmarks'][4]['Y']],
            "MidJawlineRight":[face['Landmarks'][28]['X']-face['Landmarks'][4]['X'],face['Landmarks'][28]['Y']-face['Landmarks'][4]['Y']],
            "UpperJawlineRight":[face['Landmarks'][29]['X']-face['Landmarks'][4]['X'],face['Landmarks'][29]['Y']-face['Landmarks'][4]['Y']]})
            
#create a pandas dataframe with the 'results_long' list 
img_df_long = pd.DataFrame(results_long, columns=["ID",'ChinBottom','LeftEyeDown',
'LeftEyeLeft','LeftEyeRight','LeftEyeUp','LeftEyebrowLeft','LeftEyebrowRight','LeftEyebrowUp',
'LeftPupil','MidJawlineLeft','MidJawlineRight','MouthDown','MouthLeft','MouthRight','MouthUp',
'NoseLeft','NoseRight','RightEyeDown', 'RightEyeLeft','RightEyeRight','RightEyeUp',
'RightEyebrowLeft','RightEyebrowRight','RightEyebrowUp','RightPupil','UpperJawlineLeft','UpperJawlineRight'])

#transpose it so the matrix works well with the R script
img_df_long_1 = img_df_long.transpose()

#create the csv
img_df_long_1.to_csv("schoeller_ctrl_coord_nose_norm.csv")






    

            
            
            
