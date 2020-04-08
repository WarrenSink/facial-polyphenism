for file in ~/Desktop/OneDrive-Van-Andel-Institute/Schoeller_IMGs/extracted-faces ; do
    echo mv -v "$file" "${file#*_}"
done
